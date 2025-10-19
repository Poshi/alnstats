mod bam_processor;
mod bam_stats_collector;
mod cigar_ext;
mod cli;
mod constants;
mod duplicate_stats;
mod error;
mod math;
mod statistic;
mod writer;
mod yield_stats;

use crate::bam_processor::{aggregate_stats, process_bam};
use crate::cli::Args;
use crate::error::AppError;
use crate::writer::write_results;
use clap::Parser;
use log::debug;

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