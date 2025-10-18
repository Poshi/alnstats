mod bam_stats_collector;
mod cigar_ext;
mod cli;
mod constants;
mod duplicate_stats;
mod error;
mod statistic;
mod yield_stats;

use crate::bam_stats_collector::BamStatsCollector;
use crate::cli::{Args, Aggregation};
use crate::constants::{ReadGroupTag, StatisticKind, UNKNOWN, RECORDS_LOG_INTERVAL};
use crate::error::AppError;
use clap::Parser;
use log::{debug, error, info, trace};
use noodles::bam::Record;
use noodles::bam::io::reader::Builder;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum AggregationKey {
    Sample(String),
    Library(String, String),
}

fn get_read_groups(header: &Header) -> HashMap<String, HashMap<String, String>> {
    header
        .read_groups()
        .iter()
        .map(|(k, map)| {
            let read_group_id = k.to_string().clone();

            let entry: HashMap<String, String> =
                std::iter::once((ReadGroupTag::Id.as_ref().to_string(), read_group_id.clone()))
                    .chain(
                        map.other_fields()
                            .iter()
                            .map(|(tag, value)| (format!("{tag}"), value.to_string().clone())),
                    )
                    .collect();

            (read_group_id, entry)
        })
        .collect()
}

fn get_rg_tag(record: &Record) -> Option<String> {
    record
        .data()
        .get(&Tag::READ_GROUP)
        .and_then(Result::ok)
        .and_then(|value| match value {
            Value::String(s) => Some(s.to_string()),
            _ => None, // RG should always be a String, but just in case
        })
}

type StatsPerRG = HashMap<String, BamStatsCollector>;
type StatsPerKey = HashMap<AggregationKey, BamStatsCollector>;

fn process_bam(bam_filename: &String, args: &Args) -> Result<(Header, StatsPerRG), AppError> {
    // Open input file
    trace!("Opening input file: {bam_filename}");
    let mut reader = Builder.build_from_path(bam_filename)?;

    // Read the header to position the file pointer
    trace!("Reading BAM header...");
    let header = reader.read_header()?;

    let mut stats_per_rg: StatsPerRG = HashMap::new();

    // Traverse input file while filling in the stats
    trace!("Processing BAM records...");
    for (i, rec) in reader.records().enumerate() {
        if i > 0 && i % RECORDS_LOG_INTERVAL == 0 {
            info!("{i} elements processed...");
        }

        match rec {
            Ok(record) => {
                let rg_id = get_rg_tag(&record).unwrap_or_else(|| UNKNOWN.to_string());

                let collector = stats_per_rg
                    .entry(rg_id)
                    .or_insert_with(|| BamStatsCollector::new(args));

                collector.add_record(&record);
            }
            Err(e) => {
                error!("Error reading record {i}: {e}");
            }
        }
    }

    Ok((header, stats_per_rg))
}

fn aggregate_stats(
    stats_per_rg: &StatsPerRG,
    header: &Header,
    args: &Args,
) -> HashMap<AggregationKey, BamStatsCollector> {
    let mut aggregated_stats: StatsPerKey = HashMap::new();
    let read_groups_info = get_read_groups(header);

    for (rg_id, stats_collector) in stats_per_rg {
        let rg_map = read_groups_info.get(rg_id);

        let aggregation_key = match &args.aggregation {
            Aggregation::Sample => AggregationKey::Sample(
                rg_map
                    .and_then(|rg| rg.get(ReadGroupTag::Sample.as_ref()))
                    .cloned()
                    .unwrap_or_else(|| UNKNOWN.to_string()),
            ),
            Aggregation::Library => {
                let sample = rg_map
                    .and_then(|rg| rg.get(ReadGroupTag::Sample.as_ref()))
                    .cloned()
                    .unwrap_or_else(|| UNKNOWN.to_string());
                let library = rg_map
                    .and_then(|rg| rg.get(ReadGroupTag::Library.as_ref()))
                    .cloned()
                    .unwrap_or_else(|| UNKNOWN.to_string());
                AggregationKey::Library(sample, library)
            }
        };

        let target_stats_collector = aggregated_stats
            .entry(aggregation_key)
            .or_insert_with(|| BamStatsCollector::new(args));

        *target_stats_collector += stats_collector;
    }

    aggregated_stats
}

fn write_output(
    filename: &str,
    data: &serde_json::Map<String, serde_json::Value>,
) -> Result<(), AppError> {
    trace!("Output will be written to: {filename}");
    let file = File::create(filename)?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, data)?;
    Ok(())
}

fn write_results(stats_per_aggregate_key: &StatsPerKey, args: &Args) -> Result<(), AppError> {
    trace!("Generating JSON output...");

    let mut yield_results = serde_json::Map::new();
    let mut duplicate_results = serde_json::Map::new();

    for (key, stats_collector) in stats_per_aggregate_key {
        for stat in &stats_collector.stats {
            let json_val = stat.as_json();
            match key {
                AggregationKey::Sample(sample_name) => match stat.kind() {
                    StatisticKind::YieldPE | StatisticKind::YieldSE => {
                        yield_results.insert(sample_name.clone(), json_val);
                    }
                    StatisticKind::Duplicate => {
                        duplicate_results.insert(sample_name.clone(), json_val);
                    }
                },
                AggregationKey::Library(sample_name, library_name) => {
                    let target_map = match stat.kind() {
                        StatisticKind::YieldPE | StatisticKind::YieldSE => &mut yield_results,
                        StatisticKind::Duplicate => &mut duplicate_results,
                    };

                    let sample_entry = target_map
                        .entry(sample_name.clone())
                        .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()));

                    if let Some(sample_map) = sample_entry.as_object_mut() {
                        sample_map.insert(library_name.clone(), json_val);
                    }
                }
            }
        }
    }

    if let Some(filename) = &args.yield_out {
        write_output(filename, &yield_results)?;
    }

    if let Some(filename) = &args.metrics {
        write_output(filename, &duplicate_results)?;
    }

    Ok(())
}

fn main() -> Result<(), AppError> {
    // Parse CLI arguments
    let args = Args::parse();

    // Set up logging
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();

    debug!("Arguments: {args:?}");

    // Process the BAM file, filling in the stats objects, we get stats per RGID
    let (header, stats_per_rg) = process_bam(&args.input, &args)?;

    // We aggregate those stats by sample or library
    let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);

    // Finally, we write the results to disk
    write_results(&aggregated_stats, &args)?;

    Ok(())
}
