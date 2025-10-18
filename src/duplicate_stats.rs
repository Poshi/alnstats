use log::trace;
use std::collections::HashSet;
use std::ops::AddAssign;

use noodles::bam::Record;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::value::Value;
use serde::Serialize;

use crate::constants::{DEFAULT_DUP_TAG, SEQ_DUP_VALUE, StatisticKind};
use crate::error::AppError;
use crate::statistic::Statistic;

#[derive(Debug, Clone)]
pub struct DuplicateStats {
    duplicate_type_tags: HashSet<String>,
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
            duplicate_type_tags: HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
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
            new_stats.duplicate_type_tags = HashSet::from_iter(duplicate_type_tag.iter().cloned())
        }

        new_stats
    }

    #[cfg(test)]
    pub fn unmapped_reads(&self) -> u64 {
        self.unmapped_reads
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

    /*
    pub fn total_reads(&self) -> u64 {
        (self.read_pairs_examined * 2) + self.unpaired_reads_examined + self.unmapped_reads
    }

    pub fn total_clusters(&self) -> Result<u64, RuntimeError> {
        let total_reads = self.total_reads();

        if (total_reads % 2) != 0 {
            Err(RuntimeError(String::from("The number of reads should be an even number.")))
        } else {
            Ok(total_reads / 2)
        }
    }
    */

    fn estimate_library_size(read_pairs: u64, unique_read_pairs: u64) -> Result<u64, AppError> {
        /// Estimate the size of the library.
        ///
        /// Estimates the size of a library based on the number of paired end molecules observed
        /// and the number of unique pairs observed.
        /// Based on the Lander-Waterman equation that states:
        /// C/X = 1-exp(-N/X)
        /// where
        /// X = number of distinct molecules in library
        /// N = number of read pairs
        /// C = number of distinct fragments observed in read pairs
        ///
        /// Raises:
        /// `RuntimeError`: if read pairs or duplicate pairs are zero
        ///
        /// Returns:
        /// u64: the estimated size of the library
        fn f(x: f64, c: f64, n: f64) -> f64 {
            (c / x) - 1.0 + (-n / x).exp()
        }

        fn mid(x: f64, y: f64) -> f64 {
            f64::midpoint(x, y)
        }

        let read_pair_duplicates = read_pairs - unique_read_pairs;

        if read_pairs == 0 || read_pair_duplicates == 0 {
            return Err(AppError::Runtime(String::from(
                "Read pairs or duplicate pairs are zero!",
            )));
        }

        let mut lower_bound: f64 = 1.0;
        let mut upper_bound: f64 = 100.0;

        if unique_read_pairs >= read_pairs
            || f(
                lower_bound * unique_read_pairs as f64,
                unique_read_pairs as f64,
                read_pairs as f64,
            ) < 0.0
        {
            return Err(AppError::Runtime(format!(
                "Invalid values for pairs and unique pairs: {read_pairs}, {unique_read_pairs}"
            )));
        }

        // Find value of upper_bound, large enough to act as other side for bisection method
        while f(
            upper_bound * unique_read_pairs as f64,
            unique_read_pairs as f64,
            read_pairs as f64,
        ) > 0.0
        {
            upper_bound *= 10.0;
        }

        // Use bisection method (no more than 40 times) to find solution
        for _ in 0..40 {
            let r = mid(lower_bound, upper_bound);
            let u = f(
                r * unique_read_pairs as f64,
                unique_read_pairs as f64,
                read_pairs as f64,
            );
            if u > 0.0 {
                // Move lower bound up
                lower_bound = r;
            } else if u < 0.0 {
                // Move upper bound down
                upper_bound = r;
            } else {
                break;
            }
        }

        Ok((unique_read_pairs as f64 * mid(lower_bound, upper_bound)) as u64)
    }

    pub fn estimated_library_size(&self) -> Result<u64, AppError> {
        let read_pairs = self.read_pairs_examined() - self.read_pair_optical_duplicates();
        let unique_read_pairs = self.read_pairs_examined() - self.read_pair_duplicates();

        DuplicateStats::estimate_library_size(read_pairs, unique_read_pairs)
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

    fn is_optical_duplicate(record: &Record, duplicate_type_tags: &HashSet<String>) -> bool {
        record
            .data()
            .iter()
            .filter_map(Result::ok)
            .any(|(tag, value)| {
                // Convert tag (2 bytes) to String for comparison
                let tag_str = std::str::from_utf8(tag.as_ref()).unwrap_or("");

                duplicate_type_tags.contains(tag_str)
                    && matches!(value, Value::String(s) if s == SEQ_DUP_VALUE)
            })
    }
}

#[cfg(test)]
impl DuplicateStats {
    pub(crate) fn for_test(unmapped_reads: u64) -> Self {
        Self {
            unmapped_reads,
            ..Default::default()
        }
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
            estimated_library_size: stats.estimated_library_size().unwrap_or(0),
        }
    }
}

impl Statistic for DuplicateStats {
    fn add_record(&mut self, rhs: &Record) {
        let flags = rhs.flags();

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
    }

    fn as_json(&self) -> serde_json::Value {
        serde_json::to_value(DuplicateStatsJson::from(self))
            .expect("Failed to serialize DuplicateStats to JSON")
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

impl AddAssign<&Record> for DuplicateStats {
    fn add_assign(&mut self, rhs: &Record) {
        self.add_record(rhs);
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

#[cfg(test)]
mod tests {
    use super::*;
    // use noodles::sam::alignment::record::Flags; // Removed redundant import

    #[test]
    fn test_duplicatestats_add_assign() {
        let mut stats1 = DuplicateStats {
            duplicate_type_tags: HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
            unpaired_reads_examined: 10,
            pvt_read_pairs_examined: 100,
            secondary_or_supplementary_rds: 5,
            unmapped_reads: 2,
            unpaired_read_duplicates: 1,
            pvt_read_pair_duplicates: 10,
            pvt_read_pair_optical_duplicates: 3,
        };
        let stats2 = DuplicateStats {
            duplicate_type_tags: HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
            unpaired_reads_examined: 20,
            pvt_read_pairs_examined: 200,
            secondary_or_supplementary_rds: 10,
            unmapped_reads: 4,
            unpaired_read_duplicates: 2,
            pvt_read_pair_duplicates: 20,
            pvt_read_pair_optical_duplicates: 6,
        };
        stats1 += &stats2;

        let expected_tags = HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]);

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
        assert!(stats.duplicate_type_tags.contains(DEFAULT_DUP_TAG));
        assert_eq!(stats.duplicate_type_tags.len(), 1);
    }

    #[test]
    fn test_new_custom_tag() {
        let custom_tag = String::from("XT");
        let stats = DuplicateStats::new(&[custom_tag.clone()]);
        assert!(stats.duplicate_type_tags.contains(&custom_tag));
        assert_eq!(stats.duplicate_type_tags.len(), 1);
    }

    #[test]
    fn test_new_multiple_custom_tags() {
        let tags = vec![String::from("XT"), String::from("YD")];
        let stats = DuplicateStats::new(&tags);
        assert!(stats.duplicate_type_tags.contains("XT"));
        assert!(stats.duplicate_type_tags.contains("YD"));
        assert_eq!(stats.duplicate_type_tags.len(), 2);
    }

    #[test]
    fn test_read_pairs_examined() {
        let mut stats = DuplicateStats::new(&[]);
        stats.pvt_read_pairs_examined = 100;
        assert_eq!(stats.read_pairs_examined(), 50);
    }

    #[test]
    fn test_read_pair_duplicates() {
        let mut stats = DuplicateStats::new(&[]);
        stats.pvt_read_pair_duplicates = 50;
        assert_eq!(stats.read_pair_duplicates(), 25);
    }

    #[test]
    fn test_read_pair_optical_duplicates() {
        let mut stats = DuplicateStats::new(&[]);
        stats.pvt_read_pair_optical_duplicates = 10;
        assert_eq!(stats.read_pair_optical_duplicates(), 5);
    }

    #[test]
    fn test_percent_duplication_zero_denominator() {
        let stats = DuplicateStats::new(&[]);
        assert_eq!(stats.percent_duplication(), 0.0);
    }

    #[test]
    fn test_percent_duplication_valid_cases() {
        let stats = DuplicateStats {
            duplicate_type_tags: HashSet::new(),
            unpaired_reads_examined: 10,
            pvt_read_pairs_examined: 20, // 10 pairs
            secondary_or_supplementary_rds: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 5,
            pvt_read_pair_duplicates: 4, // 2 pairs
            pvt_read_pair_optical_duplicates: 0,
        };
        // (2 * 2 + 5) / (2 * 10 + 10) = 9 / 30 = 0.3
        assert_eq!(stats.percent_duplication(), 0.3);
    }

    #[test]
    fn test_estimated_library_size_read_pairs_zero_error() {
        let stats = DuplicateStats::new(&[]);
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(err, AppError::Runtime(_)));
    }

    #[test]
    fn test_estimated_library_size_duplicate_pairs_zero_error() {
        let mut stats = DuplicateStats::new(&[]);
        stats.pvt_read_pairs_examined = 2; // 1 pair
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(err, AppError::Runtime(_)));
    }

    #[test]
    fn test_estimated_library_size_invalid_unique_pairs_error() {
        let mut stats = DuplicateStats::new(&[]);
        stats.pvt_read_pairs_examined = 2; // 1 pair
        stats.pvt_read_pair_duplicates = 0; // 0 duplicates, so unique_read_pairs = 1
        // unique_read_pairs (1) >= read_pairs (1) should trigger error
        let err = stats.estimated_library_size().unwrap_err();
        assert!(matches!(err, AppError::Runtime(_)));
    }

    #[test]
    fn test_estimated_library_size_valid_case() {
        // Example from Picard documentation
        let mut stats = DuplicateStats::new(&[]);
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
        assert!(DuplicateStats::is_secondary_or_supplementary(&Flags::SECONDARY));
        assert!(DuplicateStats::is_secondary_or_supplementary(&Flags::SUPPLEMENTARY));
        assert!(!DuplicateStats::is_secondary_or_supplementary(&Flags::empty()));
    }

    #[test]
    fn test_is_unpaired_read() {
        // Not segmented, so unpaired
        assert!(DuplicateStats::is_unpaired_read(&Flags::empty()));
        // Segmented, but mate unmapped, so unpaired
        assert!(DuplicateStats::is_unpaired_read(&(Flags::SEGMENTED | Flags::MATE_UNMAPPED)));
        // Segmented and mate mapped, so paired
        assert!(!DuplicateStats::is_unpaired_read(&(Flags::SEGMENTED)));
    }

    #[test]
    fn test_is_valid_duplicate_candidate() {
        // Is duplicate, not secondary, not supplementary, not unmapped
        assert!(DuplicateStats::is_valid_duplicate_candidate(&Flags::DUPLICATE));
        // Not duplicate
        assert!(!DuplicateStats::is_valid_duplicate_candidate(&Flags::empty()));
        // Is duplicate, but secondary
        assert!(!DuplicateStats::is_valid_duplicate_candidate(&(Flags::DUPLICATE | Flags::SECONDARY)));
        // Is duplicate, but supplementary
        assert!(!DuplicateStats::is_valid_duplicate_candidate(&(Flags::DUPLICATE | Flags::SUPPLEMENTARY)));
        // Is duplicate, but unmapped
        assert!(!DuplicateStats::is_valid_duplicate_candidate(&(Flags::DUPLICATE | Flags::UNMAPPED)));
    }

    #[test]
    fn test_duplicatestats_as_json() {
        let stats = DuplicateStats {
            duplicate_type_tags: HashSet::from_iter(vec![DEFAULT_DUP_TAG.to_string()]),
            unpaired_reads_examined: 10,
            pvt_read_pairs_examined: 100,
            secondary_or_supplementary_rds: 5,
            unmapped_reads: 2,
            unpaired_read_duplicates: 1,
            pvt_read_pair_duplicates: 10,
            pvt_read_pair_optical_duplicates: 3,
        };
        let json = stats.as_json();
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
        let stats = DuplicateStats::new(&[]);
        assert_eq!(stats.kind(), StatisticKind::Duplicate);
    }
}