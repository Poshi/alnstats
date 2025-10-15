mod statistic;
mod yield_stats;
mod duplicate_stats;
mod runtime_error;
mod cigar_ext;

use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use log::{error, info, debug, trace};
use clap::{Parser, ValueEnum};
use clap_verbosity_flag::Verbosity;
use noodles::bam::io::reader::Builder;
use noodles::bam::Record;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::yield_stats::PEYieldStats;
use crate::duplicate_stats::DuplicateStats;
use crate::statistic::Statistic;


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
        write!(f, "{}", fmt_str)
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
    #[arg(short, long="yield", value_name="YIELD")]
    yield_out: Option<String>,

    /// The tag names used for marking the duplicate type
    #[arg(short, long, default_values=&["dt"])]
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
            let read_group_id = k.to_string().to_owned();

            let entry: HashMap<String, String> = std::iter::once(("ID".to_string(), read_group_id.clone()))
                .chain(
                    map
                        .other_fields()
                        .iter()
                        .map(|(tag, value)| (format!("{}", tag), value.to_string().to_owned()))
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
        .and_then(|result| result.ok())
        .and_then(|value| match value {
            Value::String(s) => Some(s.to_string()),
            _ => None, // RG should always be a String, but just in case
        })
}

fn process_bam(bam_filename: &String, args: &Args) -> Result<(Header, HashMap<String, Vec<Box<dyn Statistic>>>), Box<dyn Error>> {
    // Open input file
    trace!("Opening input file: {}", bam_filename);
    let mut reader = Builder::default().build_from_path(bam_filename)?;

    // Read the header to position the file pointer
    trace!("Reading BAM header...");
    let header = reader.read_header()?;

    let mut stats_per_rg: HashMap<String, Vec<Box<dyn Statistic>>> = HashMap::new();

    // Traverse input file while filling in the stats
    trace!("Processing BAM records...");
    for (i, rec) in reader.records().enumerate() {
        if i > 0 && i % 10_000_000 == 0 {
            info!("{} elements processed...", i);
        }

        match rec {
            Ok(record) => {
                let rg_id = get_rg_tag(&record).unwrap_or_else(|| "unknown".to_string());

                let stats = stats_per_rg
                    .entry(rg_id)
                    .or_insert_with(|| init_stats(args));

                for stat in stats.iter_mut() {
                    stat.add_record(&record);
                }
            },
            Err(e) => {
                error!("Error reading record {}: {}", i, e);
                continue;
            }
        }
    }

    Ok((header, stats_per_rg))
}

fn aggregate_stats(stats_per_rg: &HashMap<String, Vec<Box<dyn Statistic>>>, header: &Header, args: &Args) -> HashMap<String, Vec<Box<dyn Statistic>>> {
    let mut aggregated_stats: HashMap<String, Vec<Box<dyn Statistic>>> = HashMap::new();
    let read_groups_info = get_read_groups(header);
    let aggregation_tag = match args.aggregation {
        Aggregation::Sample => "SM",
        Aggregation::Library => "LB",
    };

    for (rg_id, stats_vec) in stats_per_rg {
        let aggregation_key = read_groups_info
            .get(rg_id)
            .and_then(|rg_map| rg_map.get(aggregation_tag))
            .cloned()
            .unwrap_or_else(|| "unknown".to_string());

        let target_stats_vec = aggregated_stats
            .entry(aggregation_key)
            .or_insert_with(|| init_stats(args));

        for (i, source_stat) in stats_vec.iter().enumerate() {
            target_stats_vec[i].add_assign_to_statistic(source_stat.as_ref());
        }
    }

    aggregated_stats
}

fn write_results(stats_per_rg: &HashMap<String, Vec<Box<dyn Statistic>>>, args: &Args) -> Result<(), Box<dyn Error>> {
    trace!("Generating JSON output...");

    use std::io::BufWriter;

    let mut yield_results = serde_json::Map::new();
    let mut duplicate_results = serde_json::Map::new();

    for (rg_id, stats_vec) in stats_per_rg {
        for stat in stats_vec {
            let json_val = stat.as_json();
            match stat.kind() {
                "yield_pe" | "yield_se" => {
                    yield_results.insert(rg_id.clone(), json_val);
                }
                "duplicate" => {
                    duplicate_results.insert(rg_id.clone(), json_val);
                }
                _ => {}
            }
        }
    }

    if let Some(filename) = &args.yield_out {
        trace!("Output will be written to: {}", filename);
        let file = std::fs::File::create(filename)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &yield_results)?;
    }

    if let Some(filename) = &args.metrics {
        trace!("Output will be written to: {}", filename);
        let file = std::fs::File::create(filename)?;
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

    debug!("Arguments: {:?}", args);

    // Process the BAM file, filling in the stats objects
    match process_bam(&args.input, &args) {
        Ok((header, stats_per_rg)) => {
            let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);
            if let Err(e) = write_results(&aggregated_stats, &args) {
                error!("Error writing results: {}", e);
                std::process::exit(1);
            }
        }
        Err(e) => {
            error!("Error processing BAM file: {}", e);
            std::process::exit(1);
        }
    }
}
