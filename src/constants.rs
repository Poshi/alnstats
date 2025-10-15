// Project-wide constants.

// --- BAM Tags ---

/// Default tag for duplicate type marking.
pub const DEFAULT_DUP_TAG: &str = "dt";

/// Enum for Read Group tags.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ReadGroupTag {
    Id,
    Sample,
    Library,
}

impl AsRef<str> for ReadGroupTag {
    fn as_ref(&self) -> &str {
        match self {
            ReadGroupTag::Id => "ID",
            ReadGroupTag::Sample => "SM",
            ReadGroupTag::Library => "LB",
        }
    }
}

// --- BAM Tag Values ---

/// Value for sequencing duplicates (`duplicate_type_tag` content).
pub const SEQ_DUP_VALUE: &str = "SQ";

// --- Statistic Kinds ---

/// Enum for different kinds of statistics.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum StatisticKind {
    YieldPE,
    YieldSE,
    Duplicate,
}

impl AsRef<str> for StatisticKind {
    fn as_ref(&self) -> &str {
        match self {
            StatisticKind::YieldPE => "yield_pe",
            StatisticKind::YieldSE => "yield_se",
            StatisticKind::Duplicate => "duplicate",
        }
    }
}

// --- Default / Fallback Values ---

/// Fallback value for when a tag is not found.
pub const UNKNOWN: &str = "unknown";
