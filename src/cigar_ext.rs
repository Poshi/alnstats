use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;

pub trait CigarExt {
    fn query_alignment_length(&self) -> u32;
}

impl<C> CigarExt for C
where
    C: Cigar,
{
    fn query_alignment_length(&self) -> u32 {
        self.iter().filter_map(Result::ok).fold(0, |acc, op| {
            acc + match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion => {
                    op.len() as u32
                }
                _ => 0,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::record::Cigar;

    #[test]
    fn test_simple_match() {
        let cigar = Cigar::new(b"10M");
        assert_eq!(cigar.query_alignment_length(), 10);
    }

    #[test]
    fn test_match_and_insertion() {
        let cigar = Cigar::new(b"5M10I15M");
        assert_eq!(cigar.query_alignment_length(), 30);
    }

    #[test]
    fn test_match_and_deletion() {
        let cigar = Cigar::new(b"5M10D15M");
        assert_eq!(cigar.query_alignment_length(), 20);
    }

    #[test]
    fn test_complex_cigar() {
        let cigar = Cigar::new(b"10H20S5M1I2M4D");
        assert_eq!(cigar.query_alignment_length(), 8);
    }

    #[test]
    fn test_only_non_counting_ops() {
        let cigar = Cigar::new(b"10H20S30D40P50N");
        assert_eq!(cigar.query_alignment_length(), 0);
    }

    #[test]
    fn test_empty_cigar() {
        let cigar = Cigar::new(b"");
        assert_eq!(cigar.query_alignment_length(), 0);
    }

    #[test]
    fn test_sequence_match_and_mismatch() {
        let cigar = Cigar::new(b"5=10X15I");
        assert_eq!(cigar.query_alignment_length(), 30);
    }
}
