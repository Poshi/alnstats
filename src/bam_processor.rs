use crate::bam_stats_collector::BamStatsCollector;
use crate::cli::{Aggregation, Args};
use crate::constants::{RECORDS_LOG_INTERVAL, ReadGroupTag, UNKNOWN};
use crate::error::AppError;
use log::{info, trace, warn};
use noodles::bam::io::reader::Builder;
use noodles::bam::Record;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::Header;
use std::collections::HashMap;

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

fn get_rg_info<'a>(rg_map: Option<&'a HashMap<String, String>>, tag: &str) -> String {
    rg_map
        .and_then(|rg| rg.get(tag))
        .cloned()
        .unwrap_or_else(|| UNKNOWN.to_string())
}

pub type StatsPerRG = HashMap<String, BamStatsCollector>;
pub type StatsPerKey = HashMap<AggregationKey, BamStatsCollector>;

pub fn process_bam(bam_filename: &String, args: &Args) -> Result<(Header, StatsPerRG), AppError> {
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

                if let Err(e) = collector.add_record(&record) {
                    warn!("Error processing record: {}", e);
                }
            }
            Err(e) => {
                warn!("Error reading record {i}: {e}");
            }
        }
    }

    Ok((header, stats_per_rg))
}

pub fn aggregate_stats(
    stats_per_rg: &StatsPerRG,
    header: &Header,
    args: &Args,
) -> HashMap<AggregationKey, BamStatsCollector> {
    let mut aggregated_stats: StatsPerKey = HashMap::new();
    let read_groups_info = get_read_groups(header);

    for (rg_id, stats_collector) in stats_per_rg {
        let rg_map = read_groups_info.get(rg_id);

        let sample = get_rg_info(rg_map, ReadGroupTag::Sample.as_ref());
        let aggregation_key = match &args.aggregation {
            Aggregation::Sample => AggregationKey::Sample(sample),
            Aggregation::Library => {
                let library = get_rg_info(rg_map, ReadGroupTag::Library.as_ref());
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
