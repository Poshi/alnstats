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
