mod add_record;
mod yield_stats;
mod duplicate_stats;
mod runtime_error;

use std::fmt;
use std::fs;
use serde::Serialize;
use log::{error, info, debug, trace};
use clap::{Parser, ValueEnum};
use clap_verbosity_flag::Verbosity;
use noodles::bam::io::reader::Builder;
use noodles::bam::Record;

use crate::yield_stats::PEYieldStats;
use crate::duplicate_stats::DuplicateStats;
use crate::add_record::AddRecord;


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

#[derive(Serialize)]
enum Stats {
    Yield(PEYieldStats),
    Duplicate(DuplicateStats),
}

impl AddRecord for Stats {
    fn add_record(&mut self, record: &Record) {
        match self {
            Stats::Yield(yield_stats) => *yield_stats += record,
            Stats::Duplicate(dup_stats) => *dup_stats += record,
        };
    }
}

fn init_stats(args: &Args) -> Vec<Box<Stats>> {
    trace!("Initializing stats objects...");

    let mut stats: Vec<Box<Stats>> = Vec::new();

    match args.metrics {
        Some(_) => {
            trace!("Crating stats object for DuplicateStats");
            stats.push(Box::new(Stats::Duplicate(DuplicateStats::new(&args.duplicate_type_tag))));
        },
        None => {}
    }

    match args.yield_out {
        Some(_) => {
            trace!("Crating stats object for PEYieldStats");
            stats.push(Box::new(Stats::Yield(PEYieldStats::new())));
        },
        None => {}
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


fn process_bam(bam_filename: &String, stats: Vec<Box<Stats>>) -> Vec<Box<Stats>> {
    // Open input file
    trace!("Opening input file: {}", bam_filename);
    let mut reader = Builder::default().build_from_path(bam_filename).expect("Error reading BAM file");

    // Read the header to position the file pointer
    trace!("Reading BAM header...");
    let header = reader.read_header().expect("Issues reading header");
    let rg = get_read_groups(&header);
    debug!("Read groups found: {:?}", rg);

    // Auxiliar add_record function
    fn add_record_to_stats(mut stats: Vec<Box<Stats>>, record: &Record) -> Vec<Box<Stats>> {
        for stat in stats.iter_mut() {
            stat.add_record(record);
        }
        stats
    }

    // Traverse input file while filling in the stats
    trace!("Processing BAM records...");
    let stats = reader
        .records()
        .enumerate()
        .fold(stats, |acc, (i, x)|  {
            if i > 0 && i % 10_000_000 == 0 {
                info!("{} elements processed...", i);
            }
            add_record_to_stats(acc, &x.unwrap())
        }
    );

    stats
}

fn write_results(stats: &Vec<Box<Stats>>, args: &Args) {
    trace!("Generating JSON output...");

    for stat in stats.iter() {
        let maybe_filename = match &**stat {
            Stats::Yield(_) => &args.yield_out,
            Stats::Duplicate(_) => &args.metrics,
        };

        if let Some(filename) = maybe_filename {
            trace!("Output will be written to: {}", filename);

            let json_result = serde_json::to_string_pretty(stat).unwrap();
            fs::write(filename, json_result.as_bytes()).expect("Failed to write JSON output");
        } else {
            error!("No output file specified.");
        }
    }
}

fn main() {
    // Parse CLI arguments
    let args = Args::parse();

    // Set up logging
    env_logger::Builder::new()
        .filter_level(args.verbosity.log_level_filter())
        .init();

    debug!("Arguments: {:?}", args);

    // Create a vector with all the stats objects requested
    let stats = init_stats(&args);

    // Process the BAM file, filling in the stats objects
    let stats = process_bam(&args.input, stats);

    // Generating output for requested stats
    write_results(&stats, &args);
}
