use crate::bam_stats_collector::BamStatsCollector;
use crate::cli::{Aggregation, Args};
use crate::constants::{RECORDS_LOG_INTERVAL, ReadGroupTag, UNKNOWN};
use crate::error::AppError;
use log::{info, trace, warn};

use noodles::fasta::Repository;
use noodles::fasta::io::indexed_reader;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam::Header;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles_util::alignment::io::reader::Builder;
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

fn get_rg_tag(record: &dyn Record) -> Option<String> {
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

fn process_single_record(stats_per_rg: &mut StatsPerRG, record: &dyn Record, args: &Args) {
    let rg_id = get_rg_tag(record).unwrap_or_else(|| UNKNOWN.to_string());

    let collector = stats_per_rg
        .entry(rg_id)
        .or_insert_with(|| BamStatsCollector::new(args));

    if let Err(e) = collector.add_record(record) {
        warn!("Error processing record: {}", e);
    }
}

pub fn process_bam(aln_filename: &String, args: &Args) -> Result<(Header, StatsPerRG), AppError> {
    let mut builder = Builder::default();

    if let Some(fasta) = &args.fasta {
        let repository = indexed_reader::Builder::default()
            .build_from_path(fasta)
            .map(IndexedReader::new)
            .map(Repository::new)?;

        builder = builder.set_reference_sequence_repository(repository);
    }

    // Open input file
    trace!("Opening input file: {aln_filename}");
    let mut reader = builder.build_from_path(aln_filename)?;

    // Read the header to position the file pointer
    trace!("Reading alignment header...");
    let header = reader.read_header()?;

    let mut stats_per_rg: StatsPerRG = HashMap::new();

    // Traverse input file while filling in the stats
    trace!("Processing alignment records...");
    for (i, rec) in reader.records(&header).enumerate() {
        if i > 0 && i % RECORDS_LOG_INTERVAL == 0 {
            info!("{i} elements processed...");
        }

        match rec {
            Ok(record) => process_single_record(&mut stats_per_rg, &record, args),
            Err(e) => warn!("Error reading record {i}: {e}"),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Args;
    use crate::constants::DEFAULT_DUP_TAG;
    use crate::duplicate_stats::DuplicateStats;
    use clap::Parser;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf;
    use noodles::sam::header::record::value::map::read_group::tag;
    use noodles::sam::header::record::value::{Map, map::ReadGroup};
    use std::any::TypeId;
    use std::collections::HashSet;

    #[test]
    fn test_get_rg_tag() {
        let record = RecordBuf::builder()
            .set_data(
                [(
                    Tag::new(b'R', b'G'),
                    record_buf::data::field::Value::String("rg1".into()),
                )]
                .into_iter()
                .collect(),
            )
            .build();
        assert_eq!(get_rg_tag(&record), Some("rg1".to_string()));
    }

    #[test]
    fn test_get_rg_info() {
        let mut rg_map = HashMap::new();
        rg_map.insert("SM".to_string(), "sample1".to_string());
        assert_eq!(get_rg_info(Some(&rg_map), "SM"), "sample1".to_string());
        assert_eq!(get_rg_info(Some(&rg_map), "LB"), "unknown".to_string());
    }

    #[test]
    fn test_get_read_groups() {
        let mut header_builder = Header::builder();
        header_builder = header_builder.add_read_group("rg1", Map::<ReadGroup>::default());
        let rg2 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg2", rg2);

        let header = header_builder.build();
        let rg_info = get_read_groups(&header);

        assert_eq!(rg_info.len(), 2);
        assert_eq!(rg_info.get("rg1").unwrap().get("ID").unwrap(), "rg1");
        assert_eq!(rg_info.get("rg2").unwrap().get("SM").unwrap(), "sample1");
    }

    #[test]
    fn test_aggregate_stats_by_sample() {
        let args = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "-a",
            "sample",
        ]);
        let mut stats_per_rg = StatsPerRG::new();

        let mut collector1 = BamStatsCollector::new(&args);
        let dup_stats1 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            10,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector1
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats1));
        stats_per_rg.insert("rg1".to_string(), collector1);

        let mut collector2 = BamStatsCollector::new(&args);
        let dup_stats2 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            20,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector2
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats2));
        stats_per_rg.insert("rg2".to_string(), collector2);

        let mut header_builder = Header::builder();
        let rg1 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg1", rg1);
        let rg2 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg2", rg2);
        let header = header_builder.build();

        let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);

        assert_eq!(aggregated_stats.len(), 1);
        let key = AggregationKey::Sample("sample1".to_string());
        let collector = aggregated_stats.get(&key).unwrap();
        let dup_stats = collector
            .stats
            .get(&TypeId::of::<DuplicateStats>())
            .unwrap()
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(dup_stats.unpaired_reads_examined(), 30);
    }

    #[test]
    fn test_aggregate_stats_by_library() {
        let args = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "-a",
            "library",
        ]);
        let mut stats_per_rg = StatsPerRG::new();

        let mut collector1 = BamStatsCollector::new(&args);
        let dup_stats1 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            10,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector1
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats1));
        stats_per_rg.insert("rg1".to_string(), collector1);

        let mut collector2 = BamStatsCollector::new(&args);
        let dup_stats2 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            20,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector2
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats2));
        stats_per_rg.insert("rg2".to_string(), collector2);

        let mut header_builder = Header::builder();
        let rg1 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .insert(tag::LIBRARY, b"library1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg1", rg1);
        let rg2 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .insert(tag::LIBRARY, b"library1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg2", rg2);
        let header = header_builder.build();

        let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);

        assert_eq!(aggregated_stats.len(), 1);
        let key = AggregationKey::Library("sample1".to_string(), "library1".to_string());
        let collector = aggregated_stats.get(&key).unwrap();
        let dup_stats = collector
            .stats
            .get(&TypeId::of::<DuplicateStats>())
            .unwrap()
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(dup_stats.unpaired_reads_examined(), 30);
    }

    #[test]
    fn test_aggregate_stats_by_library_different_libraries() {
        let args = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "-a",
            "library",
        ]);
        let mut stats_per_rg = StatsPerRG::new();

        let mut collector1 = BamStatsCollector::new(&args);
        let dup_stats1 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            10,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector1
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats1));
        stats_per_rg.insert("rg1".to_string(), collector1);

        let mut collector2 = BamStatsCollector::new(&args);
        let dup_stats2 = DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.as_bytes().try_into().unwrap()]),
            20,
            0,
            0,
            0,
            0,
            0,
            0,
        );
        collector2
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats2));
        stats_per_rg.insert("rg2".to_string(), collector2);

        let mut header_builder = Header::builder();
        let rg1 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .insert(tag::LIBRARY, b"library1".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg1", rg1);
        let rg2 = Map::<ReadGroup>::builder()
            .insert(tag::SAMPLE, b"sample1".to_vec())
            .insert(tag::LIBRARY, b"library2".to_vec())
            .build()
            .unwrap();
        header_builder = header_builder.add_read_group("rg2", rg2);
        let header = header_builder.build();

        let aggregated_stats = aggregate_stats(&stats_per_rg, &header, &args);

        assert_eq!(aggregated_stats.len(), 2);
        let key = AggregationKey::Library("sample1".to_string(), "library1".to_string());
        let collector = aggregated_stats.get(&key).unwrap();
        let dup_stats = collector
            .stats
            .get(&TypeId::of::<DuplicateStats>())
            .unwrap()
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(dup_stats.unpaired_reads_examined(), 10);

        let key = AggregationKey::Library("sample1".to_string(), "library2".to_string());
        let collector = aggregated_stats.get(&key).unwrap();
        let dup_stats = collector
            .stats
            .get(&TypeId::of::<DuplicateStats>())
            .unwrap()
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(dup_stats.unpaired_reads_examined(), 20);
    }
}
