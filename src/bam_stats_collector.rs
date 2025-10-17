use crate::cli::Args;
use crate::duplicate_stats::DuplicateStats;
use crate::statistic::Statistic;
use crate::yield_stats::PEYieldStats;
use noodles::bam::Record;
use std::ops::AddAssign;

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

    pub fn process_record(&mut self, record: &Record) {
        for stat in self.stats.iter_mut() {
            stat.add_record(record);
        }
    }
}

impl AddAssign<&Self> for BamStatsCollector {
    fn add_assign(&mut self, rhs: &Self) {
        for (i, stat) in self.stats.iter_mut().enumerate() {
            stat.add_assign_to_statistic(rhs.stats[i].as_ref());
        }
    }
}
