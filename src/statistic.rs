use noodles::bam::Record;
use std::any::Any;

/// Trait to accumulate information from a BAM `Record`.
///
/// This trait is object-safe so it can be used as `Box<dyn AddRecord>`.
pub trait Statistic: Any {
    /// Accumulate a record into the stats object.
    fn add_record(&mut self, record: &Record);

    /// Return a JSON representation of the current state.
    /// Implementations should rely on `serde` to build the value.
    fn as_json(&self) -> serde_json::Value;

    /// A short static key identifying the kind of stats. Used by the caller
    /// to decide where to write the output (for example, "yield" or "duplicate").
    fn kind(&self) -> &'static str;

    /// Return self as `Any` for downcasting.
    fn as_any(&self) -> &dyn Any;

    /// Add another `AddRecord` trait object to self.
    fn add_assign_to_statistic(&mut self, other: &dyn Statistic);
}
