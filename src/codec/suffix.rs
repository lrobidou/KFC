pub use bmi::{Encoder, RevCompEncoder, RevCompIter};
use itertools::Itertools;

use crate::codec::bmi::FusedReverseIterator;

use super::bmi;

// TODO I technically could remove some allocations in this module
// By removing the need for RevCompIter

fn common_suffix_length_forward(a: u64, b: u64) -> u32 {
    // bitwise XOR: all bit to 1 indicates a difference
    let diff = a ^ b;
    // count the number of 0 before the first difference starting from the end
    let diff_pos = diff.trailing_zeros();
    // divide it by two to get the number of bases
    diff_pos / 2
}

fn common_suffix_length_for_iter<I, J>(a: I, b: J) -> usize
where
    I: Iterator<Item = u64>,
    J: Iterator<Item = u64>,
{
    let mut size_suffix = 0;
    for (me, you) in a.zip(b) {
        let suffix = common_suffix_length_forward(me, you) as usize;
        size_suffix += suffix;
        if suffix != 32 {
            return size_suffix;
        }
    }
    size_suffix
}
pub fn common_suffix_length_ff(
    forward_seq: &[u8],
    encoded_forward_seq: &[u64],
    len_encoded_seq: usize,
) -> usize {
    debug_assert!(encoded_forward_seq.len() <= (len_encoded_seq + 31) / 32);
    let encoded_self = Encoder::new(forward_seq).collect_vec();

    // I have to shift (=align) data to the right so that the prefix is easy to compute
    // TODO better way ?
    let encoded_self_aligned_rigth = FusedReverseIterator::new(&encoded_self, forward_seq.len());
    let encoded_other_aligned_rigth =
        FusedReverseIterator::new(encoded_forward_seq, len_encoded_seq);

    let size_suffix =
        common_suffix_length_for_iter(encoded_self_aligned_rigth, encoded_other_aligned_rigth);
    std::cmp::min(
        size_suffix,
        std::cmp::min(forward_seq.len(), len_encoded_seq),
    )
}

pub fn common_suffix_length_fr(
    forward_seq: &[u8],
    encoded_to_reverse: &[u64],
    len_encoded_seq: usize,
) -> usize {
    debug_assert!(encoded_to_reverse.len() <= (len_encoded_seq + 31) / 32);
    let encoded_self = Encoder::new(forward_seq).collect_vec();
    let revcomp_other = RevCompIter::new(encoded_to_reverse, len_encoded_seq).collect_vec();

    let encoded_self_aligned_rigth = FusedReverseIterator::new(&encoded_self, forward_seq.len());
    let revcomp_other_aligned_rigth = FusedReverseIterator::new(&revcomp_other, len_encoded_seq);

    let size_suffix =
        common_suffix_length_for_iter(encoded_self_aligned_rigth, revcomp_other_aligned_rigth);
    std::cmp::min(
        size_suffix,
        std::cmp::min(forward_seq.len(), len_encoded_seq),
    )
}

pub fn common_suffix_length_rf(
    seq_to_reverse: &[u8],
    encoded_forward_seq: &[u64],
    len_encoded_seq: usize,
) -> usize {
    debug_assert!(encoded_forward_seq.len() <= (len_encoded_seq + 31) / 32);
    let encoded_self = RevCompEncoder::new(seq_to_reverse).collect_vec();

    let revcomp_self_aligned_rigth = FusedReverseIterator::new(&encoded_self, seq_to_reverse.len());
    let encoded_other_aligned_rigth =
        FusedReverseIterator::new(encoded_forward_seq, len_encoded_seq);

    let size_suffix =
        common_suffix_length_for_iter(revcomp_self_aligned_rigth, encoded_other_aligned_rigth);
    std::cmp::min(
        size_suffix,
        std::cmp::min(seq_to_reverse.len(), len_encoded_seq),
    )
}

pub fn common_suffix_length_rr(
    seq_to_reverse: &[u8],
    encoded_to_reverse: &[u64],
    len_encoded_seq: usize,
) -> usize {
    debug_assert!(encoded_to_reverse.len() <= (len_encoded_seq + 31) / 32);
    let encoded_self = RevCompEncoder::new(seq_to_reverse).collect_vec();
    let revcomp_other = RevCompIter::new(encoded_to_reverse, len_encoded_seq).collect_vec();

    let revcomp_self_aligned_rigth = FusedReverseIterator::new(&encoded_self, seq_to_reverse.len());
    let encoded_other_aligned_rigth = FusedReverseIterator::new(&revcomp_other, len_encoded_seq);

    let size_suffix =
        common_suffix_length_for_iter(revcomp_self_aligned_rigth, encoded_other_aligned_rigth);
    std::cmp::min(
        size_suffix,
        std::cmp::min(seq_to_reverse.len(), len_encoded_seq),
    )
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_find_first_base_diff() {
        let a: u64 = 0b00000000_00000000_00000000_00000000_00000000_00000000_00000000_00000000;
        let b: u64 = 0b00000000_00000000_00000000_00000000_00000000_00000000_00000000_00000000;
        assert_eq!(common_suffix_length_forward(a, b), 32);
        assert_eq!(common_suffix_length_forward(b, a), 32);

        let a: u64 = 0b10000000_00000000_00000000_00000000_00000000_00000000_00000000_00000000;
        let b: u64 = 0b00000000_00000000_00000000_00000000_00000000_00000000_00000000_00000000;
        assert_eq!(common_suffix_length_forward(a, b), 31);
        assert_eq!(common_suffix_length_forward(b, a), 31);

        let a: u64 = 0b00001000_00011110_00101000_00001000_00000000_00000000_00000000_00000000;
        let b: u64 = 0b00001000_00011110_00101000_00001000_00000000_00000000_00000000_00000000;
        assert_eq!(common_suffix_length_forward(a, b), 32);
        assert_eq!(common_suffix_length_forward(b, a), 32);

        // -------------------------v------------------------v---------v-----------------------
        let a: u64 = 0b00001000_00010110_00101000_00001000_00100000_00000000_00000000_00000000;
        let b: u64 = 0b00001000_00011110_00101000_00001000_00000000_00010000_00000000_00000000;
        assert_eq!(common_suffix_length_forward(a, b), 10);
        assert_eq!(common_suffix_length_forward(b, a), 10);
    }

    #[test]
    fn test_common_suffix_length_ff() {
        {
            let a_string = String::from("ATCGGCGCATCG");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("ATCCGCGCATCG");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATGCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_fr() {
        {
            let a_string = String::from("CGATGCGCCGAT");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();
            assert_eq!(
                b,
                //  A T C G G C G C A T C G
                vec![0b00100111_11011101_00100111_00000000_00000000_00000000_00000000_00000000]
            );

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("ATCCGCGCATCG");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATGCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_rf() {
        {
            let a_string = String::from("CGATGCGCCGAT");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("CGATGCGCGGAT");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("TATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGCATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_rr() {
        {
            let a_string = String::from("CGATGCGCCGAT");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("CGATGCGCGGAT");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("TATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGCATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, b_string.len()), 0);
        }
    }
}
