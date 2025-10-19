use crate::bam_processor::{AggregationKey, StatsPerKey};
use crate::cli::Args;
use crate::error::AppError;
use crate::constants::StatisticKind;
use log::trace;
use std::fs::File;
use std::io::BufWriter;

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

pub fn write_results(stats_per_aggregate_key: &StatsPerKey, args: &Args) -> Result<(), AppError> {
    trace!("Generating JSON output...");

    let mut yield_results = serde_json::Map::new();
    let mut duplicate_results = serde_json::Map::new();

    for (key, stats_collector) in stats_per_aggregate_key {
        for stat in &stats_collector.stats {
            let json_val = stat.as_json()?;
            let target_map = match stat.kind() {
                StatisticKind::YieldPE | StatisticKind::YieldSE => &mut yield_results,
                StatisticKind::Duplicate => &mut duplicate_results,
            };
            match key {
                AggregationKey::Sample(sample_name) => {
                    target_map.insert(sample_name.clone(), json_val);
                }
                AggregationKey::Library(sample_name, library_name) => {
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
