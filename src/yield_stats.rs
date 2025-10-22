use crate::cigar_ext::CigarExt;
use crate::constants::StatisticKind;
use crate::error::AppError;
use crate::statistic::Statistic;
use log::warn;
use noodles::sam::alignment::Record;
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
    fn add_record(&mut self, record: &dyn Record) -> Result<(), AppError> {
        let seq_length = record.sequence().len() as u64;

        self.n_reads += 1;
        self.max_length = self.max_length.max(seq_length);
        self.clipped_yield += record.cigar().query_alignment_length() as u64;
        self.total_yield += seq_length;

        Ok(())
    }

    fn as_json(&self) -> Result<serde_json::Value, AppError> {
        Ok(serde_json::to_value(self)?)
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
    fn add_record(&mut self, record: &dyn Record) -> Result<(), AppError> {
        let flags = record.flags()?;

        if flags.is_supplementary() || flags.is_secondary() {
            return Ok(());
        }

        if flags.is_first_segment() {
            self.first_end.add_record(record)?;
        } else if flags.is_last_segment() {
            self.second_end.add_record(record)?;
        } else {
            warn!("Warning: read is not marked as first or last segment. Skipping.")
        }

        Ok(())
    }

    fn as_json(&self) -> Result<serde_json::Value, AppError> {
        Ok(serde_json::to_value(self)?)
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

impl AddAssign<&Self> for PEYieldStats {
    fn add_assign(&mut self, rhs: &Self) {
        self.first_end += &rhs.first_end;
        self.second_end += &rhs.second_end;
    }
}

#[cfg(test)]
impl SEYieldStats {
    pub fn n_reads(&self) -> u64 {
        self.n_reads
    }
}

#[cfg(test)]
impl PEYieldStats {
    pub fn first_end(&self) -> &SEYieldStats {
        &self.first_end
    }

    pub fn second_end(&self) -> &SEYieldStats {
        &self.second_end
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::cigar::op::Kind;

    #[test]
    fn test_seyieldstats_add_record() {
        let mut stats = SEYieldStats::default();

        let cigar = vec![Op::new(Kind::Match, 4)].into();
        let record_buf = RecordBuf::builder()
            .set_sequence(b"ACGT".to_vec().into())
            .set_cigar(cigar)
            .build();

        stats.add_record(&record_buf).unwrap();

        assert_eq!(stats.n_reads, 1);
        assert_eq!(stats.max_length, 4);
        assert_eq!(stats.clipped_yield, 4);
        assert_eq!(stats.total_yield, 4);
    }

    #[test]
    fn test_peyieldstats_add_record() {
        let mut stats = PEYieldStats::default();

        // First segment
        let cigar1 = vec![Op::new(Kind::Match, 4)].into();
        let record_buf_1 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_sequence(b"ACGT".to_vec().into())
            .set_cigar(cigar1)
            .build();
        stats.add_record(&record_buf_1).unwrap();

        // Second segment
        let cigar2 = vec![Op::new(Kind::Match, 7)].into();
        let record_buf_2 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_sequence(b"GATTACA".to_vec().into())
            .set_cigar(cigar2)
            .build();
        stats.add_record(&record_buf_2).unwrap();

        // Neither first nor last
        let record_buf_3 = RecordBuf::builder().set_flags(Flags::SEGMENTED).build();
        stats.add_record(&record_buf_3).unwrap();

        assert_eq!(stats.first_end.n_reads, 1);
        assert_eq!(stats.first_end.max_length, 4);
        assert_eq!(stats.first_end.clipped_yield, 4);
        assert_eq!(stats.first_end.total_yield, 4);

        assert_eq!(stats.second_end.n_reads, 1);
        assert_eq!(stats.second_end.max_length, 7);
        assert_eq!(stats.second_end.clipped_yield, 7);
        assert_eq!(stats.second_end.total_yield, 7);
    }

    #[test]
    fn test_seyieldstats_default() {
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
    fn test_peyieldstats_default() {
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

    #[test]
    fn test_seyieldstats_as_json() {
        let stats = SEYieldStats {
            n_reads: 1,
            max_length: 10,
            clipped_yield: 5,
            total_yield: 10,
        };
        let json = stats.as_json().unwrap();
        assert_eq!(json["n_reads"], 1);
        assert_eq!(json["max_length"], 10);
        assert_eq!(json["clipped_yield"], 5);
        assert_eq!(json["total_yield"], 10);
    }

    #[test]
    fn test_seyieldstats_kind() {
        let stats = SEYieldStats::default();
        assert_eq!(stats.kind(), StatisticKind::YieldSE);
    }

    #[test]
    fn test_peyieldstats_as_json() {
        let stats = PEYieldStats {
            first_end: SEYieldStats {
                n_reads: 1,
                max_length: 10,
                clipped_yield: 5,
                total_yield: 10,
            },
            second_end: SEYieldStats {
                n_reads: 2,
                max_length: 20,
                clipped_yield: 10,
                total_yield: 20,
            },
        };
        let json = stats.as_json().unwrap();
        assert_eq!(json["first_end"]["n_reads"], 1);
        assert_eq!(json["second_end"]["n_reads"], 2);
    }

    #[test]
    fn test_peyieldstats_kind() {
        let stats = PEYieldStats::default();
        assert_eq!(stats.kind(), StatisticKind::YieldPE);
    }

    #[test]
    fn test_seyieldstats_add_assign_to_statistic() {
        let mut stats1 = SEYieldStats {
            n_reads: 1,
            max_length: 5,
            clipped_yield: 10,
            total_yield: 15,
        };
        let stats2 = SEYieldStats {
            n_reads: 2,
            max_length: 10,
            clipped_yield: 20,
            total_yield: 30,
        };
        stats1.add_assign_to_statistic(&stats2);
        assert_eq!(stats1.n_reads, 3);
        assert_eq!(stats1.max_length, 10);
        assert_eq!(stats1.clipped_yield, 30);
        assert_eq!(stats1.total_yield, 45);
    }

    #[test]
    fn test_peyieldstats_add_assign_to_statistic() {
        let mut stats1 = PEYieldStats {
            first_end: SEYieldStats {
                n_reads: 1,
                max_length: 5,
                clipped_yield: 10,
                total_yield: 15,
            },
            second_end: SEYieldStats {
                n_reads: 1,
                max_length: 5,
                clipped_yield: 10,
                total_yield: 15,
            },
        };
        let stats2 = PEYieldStats {
            first_end: SEYieldStats {
                n_reads: 2,
                max_length: 10,
                clipped_yield: 20,
                total_yield: 30,
            },
            second_end: SEYieldStats {
                n_reads: 3,
                max_length: 12,
                clipped_yield: 22,
                total_yield: 32,
            },
        };
        stats1.add_assign_to_statistic(&stats2);
        assert_eq!(stats1.first_end.n_reads, 3);
        assert_eq!(stats1.first_end.max_length, 10);
        assert_eq!(stats1.first_end.clipped_yield, 30);
        assert_eq!(stats1.first_end.total_yield, 45);
        assert_eq!(stats1.second_end.n_reads, 4);
        assert_eq!(stats1.second_end.max_length, 12);
        assert_eq!(stats1.second_end.clipped_yield, 32);
        assert_eq!(stats1.second_end.total_yield, 47);
    }

    #[test]
    fn test_seyieldstats_add_record_with_clipping() {
        let mut stats = SEYieldStats::default();

        let cigar = vec![
            Op::new(Kind::SoftClip, 2),
            Op::new(Kind::Match, 8),
            Op::new(Kind::SoftClip, 5),
        ]
        .into();
        let record_buf = RecordBuf::builder()
            .set_sequence(b"ACGTACGTACGTACG".to_vec().into()) // 15 bases
            .set_cigar(cigar)
            .build();

        stats.add_record(&record_buf).unwrap();

        assert_eq!(stats.n_reads, 1);
        assert_eq!(stats.max_length, 15);
        // clipped_yield = 8M = 8
        assert_eq!(stats.clipped_yield, 8);
        assert_eq!(stats.total_yield, 15);
    }

    #[test]
    fn test_peyieldstats_add_record_ignore_supplementary_and_secondary() {
        let mut stats = PEYieldStats::default();

        // Supplementary record
        let record_buf_supp = RecordBuf::builder().set_flags(Flags::SUPPLEMENTARY).build();
        stats.add_record(&record_buf_supp).unwrap();

        // Secondary record
        let record_buf_sec = RecordBuf::builder().set_flags(Flags::SECONDARY).build();
        stats.add_record(&record_buf_sec).unwrap();

        assert_eq!(stats.first_end.n_reads, 0);
        assert_eq!(stats.second_end.n_reads, 0);
    }
}
