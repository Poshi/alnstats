use crate::error::AppError;

/// Auxiliar function that implements the Lander-Waterman equation
fn f(x: f64, c: f64, n: f64) -> f64 {
    (c / x) - 1.0 + (-n / x).exp()
}

/// Estimate the size of the library.
///
/// Estimates the size of a library based on the number of paired end molecules observed
/// and the number of unique pairs observed.
/// Based on the Lander-Waterman equation that states:
/// C/X = 1-exp(-N/X)
/// where
/// X = number of distinct molecules in library
/// N = number of read pairs
/// C = number of distinct fragments observed in read pairs
///
/// Raises:
/// `AppError::ZeroReads`: if read pairs or duplicate pairs are zero
/// `AppError::InvalidValues`: if there are more unique read pairs than read pairs
///
/// Returns:
/// u64: the estimated size of the library
pub fn estimate_library_size(read_pairs: u64, unique_read_pairs: u64) -> Result<u64, AppError> {
    let read_pair_duplicates = (read_pairs as i128 - unique_read_pairs as i128) as i64;

    if read_pairs == 0 || read_pair_duplicates == 0 {
        return Err(AppError::ZeroReads {
            read_pairs,
            read_pair_duplicates,
        });
    }

    let mut lower_bound: f64 = 1.0;
    let mut upper_bound: f64 = 100.0;

    if unique_read_pairs >= read_pairs
        || f(
            lower_bound * unique_read_pairs as f64,
            unique_read_pairs as f64,
            read_pairs as f64,
        ) < 0.0
    {
        return Err(AppError::InvalidValues {
            read_pairs,
            unique_read_pairs,
        });
    }

    // Find value of upper_bound, large enough to act as other side for bisection method
    while f(
        upper_bound * unique_read_pairs as f64,
        unique_read_pairs as f64,
        read_pairs as f64,
    ) > 0.0
    {
        upper_bound *= 10.0;
    }

    // Use bisection method (no more than 40 times) to find solution
    for _ in 0..40 {
        let r = f64::midpoint(lower_bound, upper_bound);
        let u = f(
            r * unique_read_pairs as f64,
            unique_read_pairs as f64,
            read_pairs as f64,
        );
        if u > 0.0 {
            // Move lower bound up
            lower_bound = r;
        } else if u < 0.0 {
            // Move upper bound down
            upper_bound = r;
        } else {
            break;
        }
    }

    Ok((unique_read_pairs as f64 * f64::midpoint(lower_bound, upper_bound)) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_matches::assert_matches;

    #[test]
    fn test_estimate_library_size_success() {
        // Typical case
        let result = estimate_library_size(1000, 800);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 2154);

        // Large values
        let result = estimate_library_size(1_000_000, 800_000);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 2154184);
    }

    #[test]
    fn test_estimate_library_size_edge_cases() {
        // unique_read_pairs very close to read_pairs
        let result = estimate_library_size(1000, 999);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 499666);
    }

    #[test]
    fn test_estimate_library_size_zero_reads() {
        // read_pairs is zero
        let result = estimate_library_size(0, 0);
        assert_matches!(
            result,
            Err(AppError::ZeroReads {
                read_pairs: 0,
                read_pair_duplicates: 0
            })
        );

        // read_pair_duplicates is zero
        let result = estimate_library_size(100, 100);
        assert_matches!(
            result,
            Err(AppError::ZeroReads {
                read_pairs: 100,
                read_pair_duplicates: 0
            })
        );
    }

    #[test]
    fn test_estimate_library_size_invalid_values() {
        // unique_read_pairs > read_pairs
        let result = estimate_library_size(100, 200);
        assert_matches!(
            result,
            Err(AppError::InvalidValues {
                read_pairs: 100,
                unique_read_pairs: 200
            })
        );

        // unique_read_pairs == read_pairs (but handled by ZeroReads)
        let result = estimate_library_size(100, 100);
        assert_matches!(
            result,
            Err(AppError::ZeroReads {
                read_pairs: 100,
                read_pair_duplicates: 0
            })
        );
    }
}
