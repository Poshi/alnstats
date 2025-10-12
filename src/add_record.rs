use noodles::bam::Record;

/// Trait to accumulate information from a BAM `Record`.
///
/// This trait is object-safe so it can be used as `Box<dyn AddRecord>`.
pub trait AddRecord {
    /// Accumulate a record into the stats object.
    fn add_record(&mut self, record: &Record);

    /// Return a JSON representation of the current state.
    /// Implementations should rely on `serde` to build the value.
    fn as_json(&self) -> serde_json::Value;

    /// A short static key identifying the kind of stats. Used by the caller
    /// to decide where to write the output (for example, "yield" or "duplicate").
    fn kind(&self) -> &'static str;
}