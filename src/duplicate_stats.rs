use std::collections::HashSet;
use std::ops::AddAssign;
use log::trace;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;
use noodles::bam::Record;
use noodles::sam::alignment::record::data::field::value::Value;
use crate::add_record::AddRecord;
use crate::runtime_error::RuntimeError;

#[derive(Debug)]
pub struct DuplicateStats {
    duplicate_type_tags: HashSet<String>,
    unpaired_reads_examined: u64,
    read_pairs_examined: u64,
    secondary_or_supplementary_rds: u64,
    unmapped_reads: u64,
    unpaired_read_duplicates: u64,
    read_pair_duplicates: u64,
    read_pair_optical_duplicates: u64,
}

impl DuplicateStats {
    pub fn new(duplicate_type_tag: &Vec<String>) -> Self {
        trace!("Creating DuplicateStats struct");

        let dt_tags: HashSet<String> = if duplicate_type_tag.is_empty() {
            HashSet::from_iter(vec![String::from("dt")])
        } else {
            HashSet::from_iter(duplicate_type_tag.iter().cloned())
        };

        DuplicateStats {
            duplicate_type_tags: dt_tags,
            unpaired_reads_examined: 0,
            read_pairs_examined: 0,
            secondary_or_supplementary_rds: 0,
            unmapped_reads: 0,
            unpaired_read_duplicates: 0,
            read_pair_duplicates: 0,
            read_pair_optical_duplicates: 0,
        }
    }

    pub fn  percent_duplication(&self) -> f64 {
        if self.read_pairs_examined == 0 || self.unpaired_reads_examined == 0 {
            return 0.0;
        }

        ((self.read_pair_duplicates * 2) + self.unpaired_read_duplicates) as f64
            / ((self.read_pairs_examined * 2) + self.unpaired_reads_examined) as f64
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

    pub fn estimated_library_size(&self) -> Result<u64, RuntimeError> {
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
    /// RuntimeError: if read pairs or duplicate pairs are zero
    ///
    /// Returns:
    /// u64: the estimated size of the library
        fn f(x: f64, c: f64, n: f64) -> f64 {
            (c / x) - 1.0 + (-n / x).exp()
        }

        fn mid(x: f64, y: f64) -> f64 {
            (x + y) / 2.0
        }

        let read_pairs = self.read_pairs_examined - self.read_pair_optical_duplicates;
        let unique_read_pairs = self.read_pairs_examined - self.read_pair_duplicates;
        let read_pair_duplicates = read_pairs - unique_read_pairs;

        if read_pairs <= 0 || read_pair_duplicates <= 0 {
            return Err(RuntimeError(String::from("Read pairs or duplicate pairs are zero!")));
        }

        let mut lower_bound: f64 = 1.0;
        let mut upper_bound: f64 = 100.0;

        if unique_read_pairs>= read_pairs || f(lower_bound * unique_read_pairs as f64, unique_read_pairs as f64, read_pairs as f64) < 0.0 {
            return Err(RuntimeError(format!("Invalid values for pairs and unique pairs: {}, {}", read_pairs, unique_read_pairs)));
        }

        // Find value of upper_bound, large enough to act as other side for bisection method
        while f(upper_bound * unique_read_pairs as f64, unique_read_pairs as f64, read_pairs as f64) > 0.0 {
            upper_bound *= 10.0;
        }

        // Use bisection method (no more than 40 times) to find solution
        for _ in 0..40 {
            let r = mid(lower_bound, upper_bound);
            let u = f(r * unique_read_pairs as f64, unique_read_pairs as f64, read_pairs as f64);
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

    pub fn to_json_with_extra(&self) -> serde_json::Value {
        trace!("Generating enriched JSON for DuplicateStats");
        serde_json::json!(
            {
                "UNPAIRED_READS_EXAMINED": self.unpaired_reads_examined,
                "READ_PAIRS_EXAMINED": self.read_pairs_examined,
                "SECONDARY_OR_SUPPLEMENTARY_RDS": self.secondary_or_supplementary_rds,
                "UNMAPPED_READS": self.unmapped_reads,
                "UNPAIRED_READ_DUPLICATES": self.unpaired_read_duplicates,
                "READ_PAIR_DUPLICATES": self.read_pair_duplicates,
                "READ_PAIR_OPTICAL_DUPLICATES": self.read_pair_optical_duplicates,
                "PERCENT_DUPLICATION": self.percent_duplication(),
                "ESTIMATED_LIBRARY_SIZE": self.estimated_library_size().unwrap_or(0),
            }
        )
    }
}

impl AddRecord for DuplicateStats {
    fn add_record(&mut self, rhs: &Record) {
        let flags = rhs.flags();

        if flags.is_supplementary() || flags.is_secondary() {
            self.secondary_or_supplementary_rds += 1;
        }

        if flags.is_unmapped() {
            self.unmapped_reads += 1;
        }

        if flags.is_mate_unmapped() {
            self.unpaired_reads_examined += 1;
            self.unpaired_read_duplicates += flags.is_duplicate() as u64;
        }

        if flags.is_first_segment() {
            self.read_pairs_examined += 1;
            if flags.is_duplicate() {
                let sq_dup = rhs
                    .data()
                    .iter()
                    .filter_map(Result::ok)
                    .any(|(tag, value)| {
                        // Convertir tag (2 bytes) a String para comparación
                        let tag_str = std::str::from_utf8(tag.as_ref()).unwrap_or("");

                        self.duplicate_type_tags.contains(tag_str) && matches!(value, Value::String(s) if s == "SQ")
                    });

                self.read_pair_duplicates += 1;
                self.read_pair_optical_duplicates += sq_dup as u64;
            }
        }
    }

    fn as_json(&self) -> serde_json::Value {
        self.to_json_with_extra()
    }

    fn kind(&self) -> &'static str { "duplicate" }
}

impl AddAssign<&Record> for DuplicateStats {
    fn add_assign(&mut self, rhs: &Record) {
        self.add_record(rhs);
    }
}

impl Serialize for DuplicateStats {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        trace!("Serializing a DuplicateStats struct");

        let mut state = serializer.serialize_struct("DuplicateStats", 1)?;
        // Here we substitute the standard serialization for its enriched JSON representation
        state.serialize_field("duplicate_stats", &self.to_json_with_extra())?;
        state.end()
    }
}
