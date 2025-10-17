use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum AppError {
    #[error("I/O error")]
    Io(#[from] io::Error),

    #[error("JSON serialization error")]
    Json(#[from] serde_json::Error),

    #[error("{0}")]
    Runtime(String),
}
