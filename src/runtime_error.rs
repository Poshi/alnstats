use std::fmt;

use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct RuntimeError(pub String);

impl fmt::Display for RuntimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RuntimeError: {}", self.0)
    }
}

impl std::error::Error for RuntimeError {}
