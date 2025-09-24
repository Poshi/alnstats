use std::ops::AddAssign;
use log::{trace, warn};
use serde::Serialize;
use noodles::bam::Record;
use crate::add_record::AddRecord;

#[derive(Debug, Serialize)]
pub struct SEYieldStats {
    n_reads: u64,
    max_length: u64,
    clipped_yield: u64,
    total_yield: u64,
}

impl SEYieldStats {
    pub fn new() -> Self {
        trace!("Creating SEYieldStats struct");
        SEYieldStats {
            n_reads: 0,
            max_length: 0,
            clipped_yield: 0,
            total_yield: 0,
        }
    }
}

impl AddRecord for SEYieldStats {
    fn add_record(&mut self, rhs: &Record) {
        let seq_length = rhs.sequence().len() as u64;
        self.n_reads += 1;
        self.max_length = self.max_length.max(seq_length);
        //self.clipped_yield += rhs.cigar().clipped_length() as u64;
        self.clipped_yield += seq_length;
        self.total_yield += seq_length;
    }
}

impl AddAssign<&Record> for SEYieldStats {
    fn add_assign(&mut self, rhs: &Record) {
        self.add_record(rhs);
    }
}

#[derive(Debug, Serialize)]
pub struct PEYieldStats {
    first_end: SEYieldStats,
    second_end: SEYieldStats,
}

impl PEYieldStats {
    pub fn new() -> Self {
        trace!("Creating PEYieldStats struct");
        PEYieldStats {
            first_end: SEYieldStats::new(),
            second_end: SEYieldStats::new(),
        }
    }
}

impl AddRecord for PEYieldStats {

    fn add_record(&mut self, rhs: &Record) {
        let flags = rhs.flags();

        if flags.is_supplementary() || flags.is_secondary() {
            return;
        }

        if flags.is_first_segment() {
            self.first_end += rhs;
        } else if flags.is_last_segment() {
            self.second_end += rhs;
        } else {
            warn!("Warning: read is not marked as first or last segment. Skipping.");
        }
    }
}
