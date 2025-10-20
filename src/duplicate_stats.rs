use log::{trace, warn};
use std::collections::HashSet;
use std::ops::AddAssign;

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::data::field::value::Value;

use noodles::sam::alignment::record::Flags;
use serde::Serialize;

use crate::constants::{DEFAULT_DUP_TAG, SEQ_DUP_VALUE, StatisticKind};
use crate::error::AppError;
use crate::math::estimate_library_size;
use crate::statistic::Statistic;

#[derive(Debug, Clone)]
pub struct DuplicateStats {
    // The pvt_* members are intended to be private.
    // They refer to "pairs", but in reality, keep an accounting of single reads.
    // When retrieving them thru the getters, their value is divided by 2 to get
    // the real value that they should contain.
    duplicate_type_tags: HashSet<[u8; 2]>,
    unpaired_reads_examined: u64,
    pvt_read_pairs_examined: u64,
    secondary_or_supplementary_rds: u64,
    unmapped_reads: u64,
    unpaired_read_duplicates: u64,
    pvt_read_pair_duplicates: u64,
    pvt_read_pair_optical_duplicates: u64,
}

impl Default for DuplicateStats {
    fn default() -> Self {
        Self {
            duplicate_type_tags: HashSet::from_iter(vec![
                DEFAULT_DUP_TAG.as_bytes().try_into().unwrap(),
            ]),
            unpaired_reads_examined: 0,
            pvt_read_pairs_examined: 0,
            secondary_or_supplementary_rds: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            pvt_read_pair_duplicates: 0,
            pvt_read_pair_optical_duplicates: 0,
        }
    }
}

impl DuplicateStats {
    pub fn new(duplicate_type_tag: &[String]) -> Self {
        trace!("Creating DuplicateStats struct");

        let mut new_stats = Self::default();

        if !duplicate_type_tag.is_empty() {
            new_stats.duplicate_type_tags = duplicate_type_tag
                .iter()
                .filter_map(|tag| {
                    if tag.len() == 2 {
                        Some(tag.as_bytes().try_into().unwrap())
                    } else {
                        warn!("Ignoring invalid duplicate tag: {tag}");
                        None
                    }
                })
                .collect()
        }

        new_stats
    }

    pub fn read_pairs_examined(&self) -> u64 {
        self.pvt_read_pairs_examined / 2
    }

    pub fn read_pair_duplicates(&self) -> u64 {
        self.pvt_read_pair_duplicates / 2
    }

    pub fn read_pair_optical_duplicates(&self) -> u64 {
        self.pvt_read_pair_optical_duplicates / 2
    }

    pub fn percent_duplication(&self) -> f64 {
        if (self.read_pairs_examined() + self.unpaired_reads_examined) == 0 {
            return 0.0;
        }

        // Originally, read_pair_duplicates and read_pairs_examined were
        // multiplied by 2 in the original Picard code.
        // But that was taking into account that they were previously
        // divided by 2 to get read pairs from single ends.
        (self.read_pair_duplicates() * 2 + self.unpaired_read_duplicates) as f64
            / (self.read_pairs_examined() * 2 + self.unpaired_reads_examined) as f64
    }

    pub fn estimated_library_size(&self) -> Result<u64, AppError> {
        let read_pairs = self.read_pairs_examined() - self.read_pair_optical_duplicates();
        let unique_read_pairs = self.read_pairs_examined() - self.read_pair_duplicates();

        estimate_library_size(read_pairs, unique_read_pairs)
    }

    fn is_unmapped_read(flags: &Flags) -> bool {
        flags.is_unmapped()
    }

    fn is_secondary_or_supplementary(flags: &Flags) -> bool {
        flags.is_supplementary() || flags.is_secondary()
    }

    fn is_unpaired_read(flags: &Flags) -> bool {
        !flags.is_segmented() || flags.is_mate_unmapped()
    }

    fn is_valid_duplicate_candidate(flags: &Flags) -> bool {
        flags.is_duplicate()
            && !(flags.is_supplementary() || flags.is_secondary() || flags.is_unmapped())
    }

    fn is_optical_duplicate(record: &dyn Record, duplicate_type_tags: &HashSet<[u8; 2]>) -> bool {
        record
            .data()
            .iter()
            .filter_map(Result::ok)
            .any(|(tag, value)| {
                duplicate_type_tags.contains(tag.as_ref())
                    && matches!(value, Value::String(s) if s == SEQ_DUP_VALUE)
            })
    }
}

impl Statistic for DuplicateStats {
    fn add_record(&mut self, rhs: &dyn Record) -> Result<(), AppError> {
        let flags = rhs.flags()?;

        if DuplicateStats::is_unmapped_read(&flags) {
            self.unmapped_reads += 1;
        } else if DuplicateStats::is_secondary_or_supplementary(&flags) {
            self.secondary_or_supplementary_rds += 1;
        } else if DuplicateStats::is_unpaired_read(&flags) {
            self.unpaired_reads_examined += 1;
        } else {
            self.pvt_read_pairs_examined += 1;
        }

        if DuplicateStats::is_valid_duplicate_candidate(&flags) {
            if DuplicateStats::is_unpaired_read(&flags) {
                self.unpaired_read_duplicates += 1;
            } else {
                self.pvt_read_pair_duplicates += 1;
                if DuplicateStats::is_optical_duplicate(rhs, &self.duplicate_type_tags) {
                    self.pvt_read_pair_optical_duplicates += 1
                }
            }
        }

        Ok(())
    }

    fn as_json(&self) -> Result<serde_json::Value, AppError> {
        Ok(serde_json::to_value(DuplicateStatsJson::from(self))?)
    }

    fn kind(&self) -> StatisticKind {
        StatisticKind::Duplicate
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn add_assign_to_statistic(&mut self, other: &dyn Statistic) {
        if let Some(other_concrete) = other.as_any().downcast_ref::<DuplicateStats>() {
            *self += other_concrete;
        }
    }
}

impl AddAssign<&Self> for DuplicateStats {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.duplicate_type_tags, rhs.duplicate_type_tags);
        self.unpaired_reads_examined += rhs.unpaired_reads_examined;
        self.pvt_read_pairs_examined += rhs.pvt_read_pairs_examined;
        self.secondary_or_supplementary_rds += rhs.secondary_or_supplementary_rds;
        self.unmapped_reads += rhs.unmapped_reads;
        self.unpaired_read_duplicates += rhs.unpaired_read_duplicates;
        self.pvt_read_pair_duplicates += rhs.pvt_read_pair_duplicates;
        self.pvt_read_pair_optical_duplicates += rhs.pvt_read_pair_optical_duplicates;
    }
}

#[derive(Debug, Serialize)]
pub struct DuplicateStatsJson {
    #[serde(rename = "UNPAIRED_READS_EXAMINED")]
    unpaired_reads_examined: u64,

    #[serde(rename = "READ_PAIRS_EXAMINED")]
    read_pairs_examined: u64,

    #[serde(rename = "SECONDARY_OR_SUPPLEMENTARY_RDS")]
    secondary_or_supplementary_rds: u64,

    #[serde(rename = "UNMAPPED_READS")]
    unmapped_reads: u64,

    #[serde(rename = "UNPAIRED_READ_DUPLICATES")]
    unpaired_read_duplicates: u64,

    #[serde(rename = "READ_PAIR_DUPLICATES")]
    read_pair_duplicates: u64,

    #[serde(rename = "READ_PAIR_OPTICAL_DUPLICATES")]
    read_pair_optical_duplicates: u64,

    #[serde(rename = "PERCENT_DUPLICATION")]
    percent_duplication: f64,

    #[serde(rename = "ESTIMATED_LIBRARY_SIZE")]
    estimated_library_size: u64,
}

impl From<&DuplicateStats> for DuplicateStatsJson {
    fn from(stats: &DuplicateStats) -> Self {
        DuplicateStatsJson {
            unpaired_reads_examined: stats.unpaired_reads_examined,
            read_pairs_examined: stats.read_pairs_examined(),
            secondary_or_supplementary_rds: stats.secondary_or_supplementary_rds,
            unmapped_reads: stats.unmapped_reads,
            unpaired_read_duplicates: stats.unpaired_read_duplicates,
            read_pair_duplicates: stats.read_pair_duplicates(),
            read_pair_optical_duplicates: stats.read_pair_optical_duplicates(),
            percent_duplication: stats.percent_duplication(),
            // We know that the calculus do not make sense for some inputs.
            // That will cause the estimated_library_size calculus to fail.
            // In that case, we want to return a zero to the user.
            estimated_library_size: stats.estimated_library_size().unwrap_or(0),
        }
    }
}

#[cfg(test)]
impl DuplicateStats {
    pub(crate) fn for_test(
        duplicate_type_tags: HashSet<[u8; 2]>,
        unpaired_reads_examined: u64,
        pvt_read_pairs_examined: u64,
        secondary_or_supplementary_rds: u64,
        unmapped_reads: u64,
        unpaired_read_duplicates: u64,
        pvt_read_pair_duplicates: u64,
        pvt_read_pair_optical_duplicates: u64,
    ) -> Self {
        Self {
            duplicate_type_tags,
            unpaired_reads_examined,
            pvt_read_pairs_examined,
            secondary_or_supplementary_rds,
            unmapped_reads,
            unpaired_read_duplicates,
            pvt_read_pair_duplicates,
            pvt_read_pair_optical_duplicates,
        }
    }

    pub fn unmapped_reads(&self) -> u64 {
        self.unmapped_reads
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use noodles::sam::alignment::record::Flags; // Removed redundant import

    #[test]
    fn test_default() {
        let stats = DuplicateStats::default();
        assert!(stats.duplicate_type_tags.contains(DEFAULT_DUP_TAG.as_bytes()));
        assert_eq!(stats.duplicate_type_tags.len(), 1);
        assert_eq!(stats.unpaired_reads_examined, 0);
        assert_eq!(stats.pvt_read_pairs_examined, 0);
        assert_eq!(stats.secondary_or_supplementary_rds, 0);
        assert_eq!(stats.unmapped_reads, 0);
        assert_eq!(stats.unpaired_read_duplicates, 0);
        assert_eq!(stats.pvt_read_pair_duplicates, 0);
        assert_eq!(stats.pvt_read_pair_optical_duplicates, 0);
    }

    #[test]
    fn test_duplicatestats_add_assign() {
        let mut stats1 = DuplicateStats::for_test(
            HashSet::from_iter(vec![*b"XT"]),
            10,
            100,
            5,
            2,
            1,
            10,
            3,
        );
        let stats2 = DuplicateStats::for_test(
            HashSet::from_iter(vec![*b"XT"]),
            20,
            200,
            10,
            4,
            2,
            20,
            6,
        );
        stats1 += &stats2;

        let expected_tags = HashSet::from_iter(vec![*b"XT"]);

        assert_eq!(stats1.duplicate_type_tags, expected_tags);
        assert_eq!(stats1.unpaired_reads_examined, 30);
        assert_eq!(stats1.pvt_read_pairs_examined, 300);
        assert_eq!(stats1.secondary_or_supplementary_rds, 15);
        assert_eq!(stats1.unmapped_reads, 6);
        assert_eq!(stats1.unpaired_read_duplicates, 3);
        assert_eq!(stats1.pvt_read_pair_duplicates, 30);
        assert_eq!(stats1.pvt_read_pair_optical_duplicates, 9);
    }

    #[test]
    fn test_new_default_tag() {
        let stats = DuplicateStats::new(&[]);
        assert!(stats.duplicate_type_tags.contains(DEFAULT_DUP_TAG.as_bytes()));
        assert_eq!(stats.duplicate_type_tags.len(), 1);
    }

    #[test]
    fn test_new_custom_tag() {
        let custom_tag = String::from("XT");
        let stats = DuplicateStats::new(&[custom_tag.clone()]);
        assert!(stats.duplicate_type_tags.contains::<[u8; 2]>(&custom_tag.as_bytes().try_into().unwrap()));
        assert_eq!(stats.duplicate_type_tags.len(), 1);
    }

    #[test]
    fn test_new_multiple_custom_tags() {
        let tags = vec![String::from("XT"), String::from("YD")];
        let stats = DuplicateStats::new(&tags);
        assert!(stats.duplicate_type_tags.contains(b"XT"));
        assert!(stats.duplicate_type_tags.contains(b"YD"));
        assert_eq!(stats.duplicate_type_tags.len(), 2);
    }

    #[test]
    fn test_read_pairs_examined() {
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pairs_examined = 100;
        assert_eq!(stats.read_pairs_examined(), 50);
    }

    #[test]
    fn test_read_pair_duplicates() {
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pair_duplicates = 50;
        assert_eq!(stats.read_pair_duplicates(), 25);
    }

    #[test]
    fn test_read_pair_optical_duplicates() {
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pair_optical_duplicates = 10;
        assert_eq!(stats.read_pair_optical_duplicates(), 5);
    }

    #[test]
    fn test_percent_duplication_zero_denominator() {
        let stats = DuplicateStats::default();
        assert_eq!(stats.percent_duplication(), 0.0);
    }

    #[test]
    fn test_percent_duplication_valid_cases() {
        let stats = DuplicateStats::for_test(
            HashSet::new(),
            10,
            20, // 10 pairs
            0,
            0,
            5,
            4, // 2 pairs
            0,
        );
        // (2 * 2 + 5) / (2 * 10 + 10) = 9 / 30 = 0.3
        assert_eq!(stats.percent_duplication(), 0.3);
    }

    #[test]
    fn test_percent_duplication_more_cases() {
        // Only unpaired reads, some duplicates
        let stats1 = DuplicateStats::for_test(HashSet::new(), 100, 0, 0, 0, 10, 0, 0);
        assert_eq!(stats1.percent_duplication(), 0.1);

        // Only paired reads, some duplicates
        let stats2 = DuplicateStats::for_test(HashSet::new(), 0, 200, 0, 0, 0, 20, 0);
        assert_eq!(stats2.percent_duplication(), 0.1);

        // No duplicates
        let stats3 = DuplicateStats::for_test(HashSet::new(), 100, 200, 0, 0, 0, 0, 0);
        assert_eq!(stats3.percent_duplication(), 0.0);

        // All duplicates
        let stats4 = DuplicateStats::for_test(HashSet::new(), 100, 200, 0, 0, 100, 200, 0);
        assert_eq!(stats4.percent_duplication(), 1.0);
    }

    #[test]
    fn test_estimated_library_size_read_pairs_zero_error() {
        let stats = DuplicateStats::default();
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(
            err,
            AppError::ZeroReads {
                read_pairs: _,
                read_pair_duplicates: _
            }
        ));
    }

    #[test]
    fn test_estimated_library_size_duplicate_pairs_zero_error() {
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pairs_examined = 2; // 1 pair
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(
            err,
            AppError::ZeroReads {
                read_pairs: _,
                read_pair_duplicates: _
            }
        ));
    }

    #[test]
    fn test_estimated_library_size_invalid_unique_pairs_error() {
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pairs_examined = 2; // 1 pair
        stats.pvt_read_pair_duplicates = 0; // 0 duplicates, so unique_read_pairs = 1
        // unique_read_pairs (1) >= read_pairs (1) should trigger error
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(
            err,
            AppError::ZeroReads {
                read_pairs: _,
                read_pair_duplicates: _
            }
        ));
    }

    #[test]
    fn test_estimated_library_size_valid_case() {
        // Example from Picard documentation
        let mut stats = DuplicateStats::default();
        stats.pvt_read_pairs_examined = 2000000; // 1,000,000 pairs
        stats.pvt_read_pair_duplicates = 1000000; // 500,000 duplicates
        // unique_read_pairs = 1,000,000 - 500,000 = 500,000
        // read_pairs = 1,000,000
        // read_pair_duplicates = 500,000
        let estimated_size = stats.estimated_library_size().unwrap();
        // The actual value from Picard is 627500. Due to floating point and bisection, it might be slightly off.
        assert!((estimated_size as f64 - 627500.0).abs() < 0.5);
    }

    use noodles::sam::alignment::record::Flags;

    #[test]
    fn test_is_unmapped_read() {
        assert!(DuplicateStats::is_unmapped_read(&Flags::UNMAPPED));
        assert!(!DuplicateStats::is_unmapped_read(&Flags::empty()));
    }

    #[test]
    fn test_is_secondary_or_supplementary() {
        assert!(DuplicateStats::is_secondary_or_supplementary(
            &Flags::SECONDARY
        ));
        assert!(DuplicateStats::is_secondary_or_supplementary(
            &Flags::SUPPLEMENTARY
        ));
        assert!(!DuplicateStats::is_secondary_or_supplementary(
            &Flags::empty()
        ));
    }

    #[test]
    fn test_is_unpaired_read() {
        // Not segmented, so unpaired
        assert!(DuplicateStats::is_unpaired_read(&Flags::empty()));
        // Segmented, but mate unmapped, so unpaired
        assert!(DuplicateStats::is_unpaired_read(
            &(Flags::SEGMENTED | Flags::MATE_UNMAPPED)
        ));
        // Segmented and mate mapped, so paired
        assert!(!DuplicateStats::is_unpaired_read(&(Flags::SEGMENTED)));
    }

    #[test]
    fn test_is_valid_duplicate_candidate() {
        // Is duplicate, not secondary, not supplementary, not unmapped
        assert!(DuplicateStats::is_valid_duplicate_candidate(
            &Flags::DUPLICATE
        ));
        // Not duplicate
        assert!(!DuplicateStats::is_valid_duplicate_candidate(
            &Flags::empty()
        ));
        // Is duplicate, but secondary
        assert!(!DuplicateStats::is_valid_duplicate_candidate(
            &(Flags::DUPLICATE | Flags::SECONDARY)
        ));
        // Is duplicate, but supplementary
        assert!(!DuplicateStats::is_valid_duplicate_candidate(
            &(Flags::DUPLICATE | Flags::SUPPLEMENTARY)
        ));
        // Is duplicate, but unmapped
        assert!(!DuplicateStats::is_valid_duplicate_candidate(
            &(Flags::DUPLICATE | Flags::UNMAPPED)
        ));
    }

    #[test]
    fn test_duplicatestats_as_json() {
        let stats = DuplicateStats::for_test(
            HashSet::from_iter(vec![*b"XT"]),
            10,
            100,
            5,
            2,
            1,
            10,
            3,
        );
        let json = stats.as_json().unwrap();
        assert_eq!(json["UNPAIRED_READS_EXAMINED"], 10);
        assert_eq!(json["READ_PAIRS_EXAMINED"], 50);
        assert_eq!(json["SECONDARY_OR_SUPPLEMENTARY_RDS"], 5);
        assert_eq!(json["UNMAPPED_READS"], 2);
        assert_eq!(json["UNPAIRED_READ_DUPLICATES"], 1);
        assert_eq!(json["READ_PAIR_DUPLICATES"], 5);
        assert_eq!(json["READ_PAIR_OPTICAL_DUPLICATES"], 1);
        // percent_duplication = (2*5 + 1) / (2*50 + 10) = 11 / 110 = 0.1
        assert_eq!(json["PERCENT_DUPLICATION"], 0.1);
        assert_eq!(json["ESTIMATED_LIBRARY_SIZE"], 283);
    }

    #[test]
    fn test_duplicatestats_kind() {
        let stats = DuplicateStats::default();
        assert_eq!(stats.kind(), StatisticKind::Duplicate);
    }
}
