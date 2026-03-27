# alnstats

`alnstats` is a high-performance command-line tool designed to calculate yield and duplicate statistics from BAM, SAM, or CRAM alignment files. Its duplicate statistics are designed to be compatible with and resemble those produced by Picard's MarkDuplicates.

## Features

- **Yield Statistics**: Computes read counts, total bases, and clipped yield for both single-end (SE) and paired-end (PE) data.
- **Duplicate Statistics**: Provides detailed metrics on duplication rates, optical duplicates, and estimated library size.
- **Aggregation**: Supports aggregating statistics at the Sample or Library level, based on the Read Group information in the alignment header.
- **Format Support**: Handles BAM, SAM, and CRAM files (CRAM requires a FASTA reference).
- **Performance**: Real-time progress updates including processing speed and interval duration.

## Usage

```bash
alnstats [OPTIONS] --input <INPUT>
```

### Parameters and Options

| Option | Short | Long | Description |
| :--- | :--- | :--- | :--- |
| **Input** | `-i` | `--input` | **Required**. The alignment file (BAM/SAM/CRAM) to process. |
| **FASTA** | `-f` | `--fasta` | The reference FASTA file. Required for decoding CRAM files. |
| **Metrics** | `-m` | `--metrics` | Output file path for duplicate metrics (JSON format). |
| **Yield** | `-y` | `--yield` | Output file path for aggregate yield results (JSON format). |
| **Tags** | `-d` | `--duplicate-type-tag` | Tag names used for marking duplicate types (default: `dt`). Can be specified multiple times. |
| **Aggregation** | `-a` | `--aggregation` | Level of data aggregation: `sample` or `library` (default: `library`). |
| **Verbosity** | `-v` | `--verbose` | Increase logging verbosity (can be used multiple times). |
| **Help** | `-h` | `--help` | Print help information. |
| **Version** | `-V` | `--version` | Print version information. |

## Output Files

`alnstats` generates two main types of output in JSON format, depending on the options provided.

### Duplicate Metrics (`--metrics`)

This file contains statistics about duplicate reads, aggregated by the chosen level (Sample or Library).

| Field | Description |
| :--- | :--- |
| `UNPAIRED_READS_EXAMINED` | Number of mapped reads examined which belong to an unpaired read or a pair where one end is unmapped. |
| `READ_PAIRS_EXAMINED` | Number of mapped read pairs examined. |
| `SECONDARY_OR_SUPPLEMENTARY_RDS` | Number of reads marked as secondary or supplementary alignments (ignored for duplicate counting). |
| `UNMAPPED_READS` | Total number of unmapped reads encountered. |
| `UNPAIRED_READ_DUPLICATES` | Number of unpaired reads marked as duplicates. |
| `READ_PAIR_DUPLICATES` | Number of read pairs marked as duplicates. |
| `READ_PAIR_OPTICAL_DUPLICATES` | Number of read pairs marked as optical duplicates (based on the provided tags). |
| `PERCENT_DUPLICATION` | The percentage of reads that are marked as duplicates. |
| `ESTIMATED_LIBRARY_SIZE` | An estimate of the number of unique molecules in the library. |

### Yield Results (`--yield`)

This file contains yield statistics for the processed alignments.

#### Paired-End (PE) Yield

For paired-end data, metrics are provided for both `first_end` and `second_end`:

| Field | Description |
| :--- | :--- |
| `n_reads` | Total number of mapped reads. |
| `max_length` | The maximum read length observed. |
| `clipped_yield` | Total number of bases that are aligned (excluding soft/hard clips). |
| `total_yield` | Total number of bases in the reads (including clips). |

#### Single-End (SE) Yield

For single-end data, the same fields (`n_reads`, `max_length`, `clipped_yield`, `total_yield`) are provided at the root of the aggregation key.

## Requirements

- Rust 1.89.0 or newer.
- For CRAM processing, the corresponding reference FASTA file must be available.

## Installation

To build `alnstats` from source:

```bash
cargo build --release
```

The resulting binary will be located at `target/release/alnstats`.
