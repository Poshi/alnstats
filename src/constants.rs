// Project-wide constants.

// --- BAM Tags ---

/// Default tag for duplicate type marking.
pub const DEFAULT_DUP_TAG: &str = "dt";
/// Tag for Sample ID in Read Group header.
pub const RG_SAMPLE_TAG: &str = "SM";
/// Tag for Library ID in Read Group header.
pub const RG_LIBRARY_TAG: &str = "LB";
/// Tag for Read Group ID in header.
pub const RG_ID_TAG: &str = "ID";

// --- BAM Tag Values ---

/// Value for sequencing duplicates (`duplicate_type_tag` content).
pub const SEQ_DUP_VALUE: &str = "SQ";

// --- Statistic Kinds ---

/// Kind identifier for Paired-End Yield stats.
pub const KIND_YIELD_PE: &str = "yield_pe";
/// Kind identifier for Single-End Yield stats.
pub const KIND_YIELD_SE: &str = "yield_se";
/// Kind identifier for Duplicate stats.
pub const KIND_DUPLICATE: &str = "duplicate";

// --- Default / Fallback Values ---

/// Fallback value for when a tag is not found.
pub const UNKNOWN: &str = "unknown";
