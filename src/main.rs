mod add_record;
mod yield_stats;
mod duplicate_stats;
mod runtime_error;

use std::collections::HashMap;
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

fn init_stats(args: &Args) -> Vec<Box<dyn AddRecord>> {
    trace!("Initializing stats objects...");

    let mut stats: Vec<Box<dyn AddRecord>> = Vec::new();

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
            _ => None, // RG debería ser siempre String, pero por si acaso
        })
}

fn process_bam(bam_filename: &String, mut stats: Vec<Box<dyn AddRecord>>) -> Vec<Box<dyn AddRecord>> {
    // Open input file
    trace!("Opening input file: {}", bam_filename);
    let mut reader = Builder::default().build_from_path(bam_filename).expect("Error reading BAM file");

    // Read the header to position the file pointer
    trace!("Reading BAM header...");
    let header = reader.read_header().expect("Issues reading header");
    let rg = get_read_groups(&header);
    debug!("Read groups found: {:?}", rg);

    // Auxiliar add_record function
    fn add_record_to_stats(stats: &mut [Box<dyn AddRecord>], record: &Record) {
        /*
        let rg = get_rg_tag(&record);
        let rg_str = match rg {
            Some(rg_id) => rg_id,
            None => String::from(""),
        };
        trace!("RG: {rg_str}");
        */

        for stat in stats.iter_mut() {
            stat.add_record(record);
        }
    }

    // Traverse input file while filling in the stats
    trace!("Processing BAM records...");
    for (i, rec) in reader.records().enumerate() {
        if i > 0 && i % 10_000_000 == 0 {
            info!("{} elements processed...", i);
        }

        match rec {
            Ok(record) => add_record_to_stats(&mut stats, &record),
            Err(e) => {
                error!("Error reading record {}: {}", i, e);
                continue;
            }
        }
    }

    stats
}

fn write_results(stats: &[Box<dyn AddRecord>], args: &Args) {
    trace!("Generating JSON output...");

    use std::io::BufWriter;

    for stat in stats.iter() {
        let kind = stat.kind();
        let maybe_filename = match kind {
            "yield_pe" | "yield_se" => &args.yield_out,
            "duplicate" => &args.metrics,
            _ => &None,
        };

        if let Some(filename) = maybe_filename {
            trace!("Output will be written to: {}", filename);
            let file = std::fs::File::create(filename).expect("Failed to create output file");
            let writer = BufWriter::new(file);
            let json_value = stat.as_json();
            serde_json::to_writer_pretty(writer, &json_value).expect("Failed to write JSON output");
        } else {
            error!("No output file specified for kind {}.", kind);
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
