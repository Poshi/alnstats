use crate::constants::DEFAULT_DUP_TAG;
use clap::{Parser, ValueEnum};
use std::fmt;

#[derive(Debug, Clone, ValueEnum)]
pub enum Aggregation {
    Sample,
    Library,
}

impl fmt::Display for Aggregation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let fmt_str = match self {
            Aggregation::Sample => "sample",
            Aggregation::Library => "library",
        };
        write!(f, "{fmt_str}")
    }
}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum AggregationKey {
    Sample(String),
    Library(String, String),
}

/// Simple program to greet a person
#[derive(Debug, Parser)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// The alignment file to process
    #[arg(short, long)]
    pub input: String,

    /// Output file for the duplicate metrics
    #[arg(short, long)]
    pub metrics: Option<String>,

    /// Output file for the aggregate yield results
    #[arg(short, long = "yield", value_name = "YIELD")]
    pub yield_out: Option<String>,

    /// The tag names used for marking the duplicate type
    #[arg(short, long, default_values=&[DEFAULT_DUP_TAG])]
    pub duplicate_type_tag: Vec<String>,

    /// Level at which we want the data to be aggregated
    #[arg(short, long, default_value_t=Aggregation::Library)]
    pub aggregation: Aggregation,

    /// Be verbose
    #[command(flatten)]
    pub verbosity: clap_verbosity_flag::Verbosity,
}
