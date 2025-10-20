use crate::error::AppError;
use std::any::TypeId;
use std::collections::HashMap;
use std::ops::AddAssign;

use noodles::sam::alignment::Record;

use crate::cli::Args;
use crate::duplicate_stats::DuplicateStats;
use crate::statistic::Statistic;
use crate::yield_stats::PEYieldStats;

pub struct BamStatsCollector {
    pub stats: HashMap<TypeId, Box<dyn Statistic>>,
}

impl BamStatsCollector {
    pub fn new(args: &Args) -> Self {
        let mut stats: HashMap<TypeId, Box<dyn Statistic>> = HashMap::new();

        if args.metrics.is_some() {
            stats.insert(
                TypeId::of::<DuplicateStats>(),
                Box::new(DuplicateStats::new(&args.duplicate_type_tag)),
            );
        }

        if args.yield_out.is_some() {
            stats.insert(
                TypeId::of::<PEYieldStats>(),
                Box::new(PEYieldStats::default()),
            );
        }

        Self { stats }
    }

    pub fn add_record(&mut self, record: &dyn Record) -> Result<(), AppError> {
        for stat in self.stats.values_mut() {
            stat.add_record(record)?
        }

        Ok(())
    }
}

impl AddAssign<&Self> for BamStatsCollector {
    fn add_assign(&mut self, rhs: &Self) {
        for (type_id, stat) in &rhs.stats {
            if let Some(self_stat) = self.stats.get_mut(type_id) {
                self_stat.add_assign_to_statistic(stat.as_ref());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Args;

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
        assert!(
            collector
                .stats
                .contains_key(&TypeId::of::<DuplicateStats>())
        );
    }

    #[test]
    fn test_new_with_yield() {
        let args = Args::parse_from(&["bamstats", "-i", "test.bam", "--yield", "yield.txt"]);
        let collector = BamStatsCollector::new(&args);
        assert_eq!(collector.stats.len(), 1);
        assert!(collector.stats.contains_key(&TypeId::of::<PEYieldStats>()));
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
        assert!(
            collector
                .stats
                .contains_key(&TypeId::of::<DuplicateStats>())
        );
        assert!(collector.stats.contains_key(&TypeId::of::<PEYieldStats>()));
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

        let dup_stats1 =
            DuplicateStats::for_test(HashSet::from_iter(vec![*b"XT"]), 0, 0, 0, 10, 0, 0, 0);
        collector1
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats1));

        let dup_stats2 =
            DuplicateStats::for_test(HashSet::from_iter(vec![*b"XT"]), 0, 0, 0, 20, 0, 0, 0);
        collector2
            .stats
            .insert(TypeId::of::<DuplicateStats>(), Box::new(dup_stats2));

        collector1 += &collector2;

        let final_stats = collector1.stats[&TypeId::of::<DuplicateStats>()]
            .as_any()
            .downcast_ref::<DuplicateStats>()
            .unwrap();
        assert_eq!(final_stats.unmapped_reads(), 30);
    }

    #[test]
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
