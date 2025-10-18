// Project-wide constants.

use strum::{AsRefStr, Display};

// --- BAM Tags ---

/// Default tag for duplicate type marking.
pub const DEFAULT_DUP_TAG: &str = "dt";

/// Enum for Read Group tags.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Display, AsRefStr)]
pub enum ReadGroupTag {
    #[strum(serialize = "ID")]
    Id,
    #[strum(serialize = "SM")]
    Sample,
    #[strum(serialize = "LB")]
    Library,
}

// --- BAM Tag Values ---

/// Value for sequencing duplicates (`duplicate_type_tag` content).
pub const SEQ_DUP_VALUE: &str = "SQ";

// --- Statistic Kinds ---

/// Enum for different kinds of statistics.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Display, AsRefStr)]
#[strum(serialize_all = "snake_case")]
pub enum StatisticKind {
    YieldPE,
    YieldSE,
    Duplicate,
}


// --- Default / Fallback Values ---

/// Fallback value for when a tag is not found.
pub const UNKNOWN: &str = "unknown";
