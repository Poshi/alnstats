mod cigar_ext;
mod constants;
mod duplicate_stats;
mod runtime_error;
mod statistic;
mod yield_stats;

use clap::{Parser, ValueEnum};
use clap_verbosity_flag::Verbosity;
use log::{debug, error, info, trace};
use noodles::bam::Record;
use noodles::bam::io::reader::Builder;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::BufWriter;

use crate::constants::{
    DEFAULT_DUP_TAG, KIND_DUPLICATE, KIND_YIELD_PE, KIND_YIELD_SE, RG_ID_TAG, RG_LIBRARY_TAG,
    RG_SAMPLE_TAG, UNKNOWN,
};
use crate::duplicate_stats::DuplicateStats;
use crate::statistic::Statistic;
use crate::yield_stats::PEYieldStats;

#[derive(Debug, Clone, ValueEnum)]
enum Aggregation {
    Sample,
    Library,
}

impl fmt::Display for Aggregation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let fmt_str = match self {
            Aggregation::Sample => "sample",
            Aggregation::Library => "library",
        };
        write!(f, "{fmt_str}")
    }
}

/// Simple program to greet a person
#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
struct Args {
    /// The alignment file to process
    #[arg(short, long)]
    input: String,

    /// Output file for the duplicate metrics
    #[arg(short, long)]
    metrics: Option<String>,

    /// Output file for the aggregate yield results
    #[arg(short, long = "yield", value_name = "YIELD")]
    yield_out: Option<String>,

    /// The tag names used for marking the duplicate type
    #[arg(short, long, default_values=&[DEFAULT_DUP_TAG])]
    duplicate_type_tag: Vec<String>,

    /// Level at which we want the data to be aggregated
    #[arg(short, long, default_value_t=Aggregation::Library)]
    aggregation: Aggregation,

    /// Be verbose
    #[command(flatten)]
    verbosity: Verbosity,
}

fn init_stats(args: &Args) -> Vec<Box<dyn Statistic>> {
    trace!("Initializing stats objects...");

    let mut stats: Vec<Box<dyn Statistic>> = Vec::new();

    if args.metrics.is_some() {
        trace!("Creating stats object for DuplicateStats");
        stats.push(Box::new(DuplicateStats::new(&args.duplicate_type_tag)));
    }

    if args.yield_out.is_some() {
        trace!("Creating stats object for PEYieldStats");
        stats.push(Box::new(PEYieldStats::new()));
    }

    stats
}

fn get_read_groups(header: &Header) -> HashMap<String, HashMap<String, String>> {
    header
        .read_groups()
        .iter()
        .map(|(k, map)| {
            let read_group_id = k.to_string().clone();

            let entry: HashMap<String, String> =
                std::iter::once((RG_ID_TAG.to_string(), read_group_id.clone()))
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

type StatsPerKey = HashMap<String, Vec<Box<dyn Statistic>>>;

fn process_bam(
    bam_filename: &String,
    args: &Args,
) -> Result<(Header, StatsPerKey), Box<dyn Error>> {
    // Open input file
    trace!("Opening input file: {bam_filename}");
    let mut reader = Builder.build_from_path(bam_filename)?;

    // Read the header to position the file pointer
    trace!("Reading BAM header...");
    let header = reader.read_header()?;

    let mut stats_per_rg: StatsPerKey = HashMap::new();

    // Traverse input file while filling in the stats
    trace!("Processing BAM records...");
    for (i, rec) in reader.records().enumerate() {
        if i > 0 && i % 10_000_000 == 0 {
            info!("{i} elements processed...");
        }

        match rec {
            Ok(record) => {
                let rg_id = get_rg_tag(&record).unwrap_or_else(|| UNKNOWN.to_string());

                let stats = stats_per_rg
                    .entry(rg_id)
                    .or_insert_with(|| init_stats(args));

                for stat in stats.iter_mut() {
                    stat.add_record(&record);
                }
            }
            Err(e) => {
                error!("Error reading record {i}: {e}");
            }
        }
    }

    Ok((header, stats_per_rg))
}

fn aggregate_stats(stats_per_rg: &StatsPerKey, header: &Header, args: &Args) -> StatsPerKey {
    let mut aggregated_stats: StatsPerKey = HashMap::new();
    let read_groups_info = get_read_groups(header);

    for (rg_id, stats_vec) in stats_per_rg {
        let rg_map = read_groups_info.get(rg_id);

        let aggregation_key = match args.aggregation {
            Aggregation::Sample => rg_map
                .and_then(|rg| rg.get(RG_SAMPLE_TAG))
                .cloned()
                .unwrap_or_else(|| UNKNOWN.to_string()),
            Aggregation::Library => {
                let sample = rg_map
                    .and_then(|rg| rg.get(RG_SAMPLE_TAG))
                    .cloned()
                    .unwrap_or_else(|| UNKNOWN.to_string());
                let library = rg_map
                    .and_then(|rg| rg.get(RG_LIBRARY_TAG))
                    .cloned()
                    .unwrap_or_else(|| UNKNOWN.to_string());
                format!("{sample} {library}")
            }
        };

        let target_stats_vec = aggregated_stats
            .entry(aggregation_key)
            .or_insert_with(|| init_stats(args));

        for (i, source_stat) in stats_vec.iter().enumerate() {
            target_stats_vec[i].add_assign_to_statistic(source_stat.as_ref());
        }
    }

    aggregated_stats
}

fn write_results(stats_per_aggregate_key: &StatsPerKey, args: &Args) -> Result<(), Box<dyn Error>> {
    trace!("Generating JSON output...");

    let mut yield_results = serde_json::Map::new();
    let mut duplicate_results = serde_json::Map::new();

    for (key, stats_vec) in stats_per_aggregate_key {
        for stat in stats_vec {
            let json_val = stat.as_json();
            match args.aggregation {
                Aggregation::Sample => match stat.kind() {
                    KIND_YIELD_PE | KIND_YIELD_SE => {
                        yield_results.insert(key.clone(), json_val);
                    }
                    KIND_DUPLICATE => {
                        duplicate_results.insert(key.clone(), json_val);
                    }
                    _ => {}
                },
                Aggregation::Library => {
                    let parts: Vec<&str> = key.splitn(2, ' ').collect();
                    let sample_name = parts[0].to_string();
                    let library_name = parts[1].to_string();

                    let target_map = match stat.kind() {
                        KIND_YIELD_PE | KIND_YIELD_SE => &mut yield_results,
                        KIND_DUPLICATE => &mut duplicate_results,
                        _ => continue,
                    };

                    let sample_entry = target_map
                        .entry(sample_name)
                        .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()));

                    if let Some(sample_map) = sample_entry.as_object_mut() {
                        sample_map.insert(library_name, json_val);
                    }
                }
            }
        }
    }

    if let Some(filename) = &args.yield_out {
        trace!("Output will be written to: {filename}");
        let file = File::create(filename)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &yield_results)?;
    }

    if let Some(filename) = &args.metrics {
        trace!("Output will be written to: {filename}");
        let file = File::create(filename)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &duplicate_results)?;
    }

    Ok(())
}

fn main() {
    // Parse CLI arguments
    let args = Args::parse();

    // Set up logging
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();

    debug!("Arguments: {args:?}");

    // Process the BAM file, filling in the stats objects
    match process_bam(&args.input, &args) {
        Ok((header, stats_per_rg)) => {
            let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);
            if let Err(e) = write_results(&aggregated_stats, &args) {
                error!("Error writing results: {e}");
                std::process::exit(1);
            }
        }
        Err(e) => {
            error!("Error processing BAM file: {e}");
            std::process::exit(1);
        }
    }
}
