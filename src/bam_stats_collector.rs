use crate::error::AppError;
use std::ops::AddAssign;

use log::warn;
use noodles::bam::Record;

use crate::cli::Args;
use crate::duplicate_stats::DuplicateStats;
use crate::statistic::Statistic;
use crate::yield_stats::PEYieldStats;

pub struct BamStatsCollector {
    pub stats: Vec<Box<dyn Statistic>>,
}

impl BamStatsCollector {
    pub fn new(args: &Args) -> Self {
        let mut stats: Vec<Box<dyn Statistic>> = Vec::new();

        if args.metrics.is_some() {
            stats.push(Box::new(DuplicateStats::new(&args.duplicate_type_tag)));
        }

        if args.yield_out.is_some() {
            stats.push(Box::new(PEYieldStats::default()));
        }

        Self { stats }
    }

    pub fn add_record(&mut self, record: &Record) -> Result<(), AppError> {
        for stat in self.stats.iter_mut() {
            if let Err(error) = stat.add_record(record) {
                match error {
                    AppError::NotFirstNotLastSegment() => {
                        warn!("Warning: read is not marked as first or last segment. Skipping.");
                        return Ok(())
                    }
                    _ => return Err(error)
                }
            }
        }

        Ok(())
    }
}

impl AddAssign<&Self> for BamStatsCollector {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.stats.len(),
            rhs.stats.len(),
            "Cannot merge BamStatsCollectors with different sets of statistics."
        );
        for (i, stat) in self.stats.iter_mut().enumerate() {
            stat.add_assign_to_statistic(rhs.stats[i].as_ref());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Args;
    use crate::constants::DEFAULT_DUP_TAG;
    use clap::Parser;
    use std::collections::HashSet;

    #[test]
    fn test_new_no_stats() {
        let args = Args::parse_from(&["bamstats", "-i", "test.bam"]);
        let collector = BamStatsCollector::new(&args);
        assert!(collector.stats.is_empty());
    }

    #[test]
    fn test_new_with_metrics() {
        let args = Args::parse_from(&["bamstats", "-i", "test.bam", "-m", "metrics.txt"]);
        let collector = BamStatsCollector::new(&args);
        assert_eq!(collector.stats.len(), 1);
        assert!(collector.stats[0]
            .as_any()
            .is::<DuplicateStats>());
    }

    #[test]
    fn test_new_with_yield() {
        let args = Args::parse_from(&["bamstats", "-i", "test.bam", "--yield", "yield.txt"]);
        let collector = BamStatsCollector::new(&args);
        assert_eq!(collector.stats.len(), 1);
        assert!(collector.stats[0].as_any().is::<PEYieldStats>());
    }

    #[test]
    fn test_new_with_all_stats() {
        let args = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "--yield",
            "yield.txt",
        ]);
        let collector = BamStatsCollector::new(&args);
        assert_eq!(collector.stats.len(), 2);
        assert!(collector.stats[0]
            .as_any()
            .is::<DuplicateStats>());
        assert!(collector.stats[1].as_any().is::<PEYieldStats>());
    }

    #[test]
    fn test_add_assign() {
        let args = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "--yield",
            "yield.txt",
        ]);
        let mut collector1 = BamStatsCollector::new(&args);
        let mut collector2 = BamStatsCollector::new(&args);

        collector1.stats[0] = Box::new(DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
            0,
            0,
            0,
            10,
            0,
            0,
            0,
        ));
        collector2.stats[0] = Box::new(DuplicateStats::for_test(
            HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
            0,
            0,
            0,
            20,
            0,
            0,
            0,
        ));

        collector1 += &collector2;

        let final_stats = collector1.stats[0]
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(final_stats.unmapped_reads(), 30);
    }

    #[test]
    #[should_panic]
    fn test_add_assign_panic() {
        let args1 = Args::parse_from(&["bamstats", "-i", "test.bam", "-m", "metrics.txt"]);
        let mut collector1 = BamStatsCollector::new(&args1);

        let args2 = Args::parse_from(&[
            "bamstats",
            "-i",
            "test.bam",
            "-m",
            "metrics.txt",
            "--yield",
            "yield.txt",
        ]);
        let collector2 = BamStatsCollector::new(&args2);

        collector1 += &collector2;
    }
}