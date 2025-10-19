use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("I/O error")]
    Io(#[from] io::Error),

    #[error("JSON serialization error")]
    Json(#[from] serde_json::Error),

    #[error("Read neither first nor last in template")]
    NotFirstNotLastSegment(),

    #[error("Read pairs ({read_pairs}) or duplicate pairs ({read_pair_duplicates}) are zero!")]
    ZeroReads {
        read_pairs: u64,
        read_pair_duplicates: u64,
    },

    #[error("Invalid values for pairs and unique pairs: {read_pairs}, {unique_read_pairs}")]
    InvalidValues {
        read_pairs: u64,
        unique_read_pairs: u64,
    },
}
