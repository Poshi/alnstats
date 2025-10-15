use crate::cigar_ext::CigarExt;
use crate::constants::StatisticKind;
use crate::statistic::Statistic;
use log::warn;
use noodles::bam::Record;
use serde::Serialize;
use std::ops::AddAssign;

#[derive(Debug, Serialize, PartialEq, Clone, Default)]
pub struct SEYieldStats {
    n_reads: u64,
    max_length: u64,
    clipped_yield: u64,
    total_yield: u64,
}

impl Statistic for SEYieldStats {
    fn add_record(&mut self, record: &Record) {
        let seq_length = record.sequence().len() as u64;

        self.n_reads += 1;
        self.max_length = self.max_length.max(seq_length);
        self.clipped_yield += u64::from(record.cigar().query_alignment_length());
        self.total_yield += seq_length;
    }

    fn as_json(&self) -> serde_json::Value {
        serde_json::to_value(self).expect("Failed to serialize SEYieldStats to JSON")
    }

    fn kind(&self) -> StatisticKind {
        StatisticKind::YieldSE
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn add_assign_to_statistic(&mut self, other: &dyn Statistic) {
        if let Some(other_concrete) = other.as_any().downcast_ref::<SEYieldStats>() {
            *self += other_concrete;
        }
    }
}

impl AddAssign<&Record> for SEYieldStats {
    fn add_assign(&mut self, rhs: &Record) {
        self.add_record(rhs);
    }
}

impl AddAssign<&Self> for SEYieldStats {
    fn add_assign(&mut self, rhs: &Self) {
        self.n_reads += rhs.n_reads;
        self.max_length = self.max_length.max(rhs.max_length);
        self.clipped_yield += rhs.clipped_yield;
        self.total_yield += rhs.total_yield;
    }
}

#[derive(Debug, Serialize, PartialEq, Clone, Default)]
pub struct PEYieldStats {
    first_end: SEYieldStats,
    second_end: SEYieldStats,
}

impl Statistic for PEYieldStats {
    fn add_record(&mut self, record: &Record) {
        let flags = record.flags();

        if flags.is_supplementary() || flags.is_secondary() {
            return;
        }

        if flags.is_first_segment() {
            self.first_end += record;
        } else if flags.is_last_segment() {
            self.second_end += record;
        } else {
            warn!("Warning: read is not marked as first or last segment. Skipping.");
        }
    }

    fn as_json(&self) -> serde_json::Value {
        serde_json::to_value(self).expect("Failed to serialize PEYieldStats to JSON")
    }

    fn kind(&self) -> StatisticKind {
        StatisticKind::YieldPE
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn add_assign_to_statistic(&mut self, other: &dyn Statistic) {
        if let Some(other_concrete) = other.as_any().downcast_ref::<PEYieldStats>() {
            *self += other_concrete;
        }
    }
}

impl AddAssign<&Record> for PEYieldStats {
    fn add_assign(&mut self, rhs: &Record) {
        self.add_record(rhs);
    }
}

impl AddAssign<&Self> for PEYieldStats {
    fn add_assign(&mut self, rhs: &Self) {
        self.first_end += &rhs.first_end;
        self.second_end += &rhs.second_end;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_seyieldstats() {
        let result = SEYieldStats::default();
        assert_eq!(
            result,
            SEYieldStats {
                n_reads: 0,
                max_length: 0,
                clipped_yield: 0,
                total_yield: 0,
            }
        );
    }

    #[test]
    fn new_peyieldstats() {
        let result = PEYieldStats::default();
        assert_eq!(
            result,
            PEYieldStats {
                first_end: SEYieldStats {
                    n_reads: 0,
                    max_length: 0,
                    clipped_yield: 0,
                    total_yield: 0,
                },
                second_end: SEYieldStats {
                    n_reads: 0,
                    max_length: 0,
                    clipped_yield: 0,
                    total_yield: 0,
                },
            }
        );
    }

    #[test]
    fn test_seyieldstats_add_assign() {
        let mut stats1 = SEYieldStats {
            n_reads: 10,
            max_length: 100,
            clipped_yield: 5,
            total_yield: 1000,
        };
        let stats2 = SEYieldStats {
            n_reads: 20,
            max_length: 150,
            clipped_yield: 10,
            total_yield: 2000,
        };
        stats1 += &stats2;
        assert_eq!(
            stats1,
            SEYieldStats {
                n_reads: 30,
                max_length: 150,
                clipped_yield: 15,
                total_yield: 3000,
            }
        );
    }

    #[test]
    fn test_peyieldstats_add_assign() {
        let mut stats1 = PEYieldStats {
            first_end: SEYieldStats {
                n_reads: 10,
                max_length: 100,
                clipped_yield: 5,
                total_yield: 1000,
            },
            second_end: SEYieldStats {
                n_reads: 10,
                max_length: 100,
                clipped_yield: 5,
                total_yield: 1000,
            },
        };
        let stats2 = PEYieldStats {
            first_end: SEYieldStats {
                n_reads: 20,
                max_length: 150,
                clipped_yield: 10,
                total_yield: 2000,
            },
            second_end: SEYieldStats {
                n_reads: 20,
                max_length: 150,
                clipped_yield: 10,
                total_yield: 2000,
            },
        };
        stats1 += &stats2;
        assert_eq!(
            stats1,
            PEYieldStats {
                first_end: SEYieldStats {
                    n_reads: 30,
                    max_length: 150,
                    clipped_yield: 15,
                    total_yield: 3000,
                },
                second_end: SEYieldStats {
                    n_reads: 30,
                    max_length: 150,
                    clipped_yield: 15,
                    total_yield: 3000,
                },
            }
        );
    }
}
