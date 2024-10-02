use crate::codec::bmi::{
    EncoderFromTheEndRightAligned, FusedReverseIterator, RevCompEncoderSRA, RevCompIterSRA,
};

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
    start: usize,
    end: usize,
) -> usize {
    let len_encoded = end - start;

    debug_assert!(encoded_forward_seq.len() >= (len_encoded + 31) / 32);
    // I have to shift (=align) data to the right so that the prefix is easy to compute
    let encoded_self_aligned_rigth = EncoderFromTheEndRightAligned::new(forward_seq);
    let encoded_other_aligned_rigth = FusedReverseIterator::new(encoded_forward_seq, start, end);

    let size_suffix =
        common_suffix_length_for_iter(encoded_self_aligned_rigth, encoded_other_aligned_rigth);
    std::cmp::min(size_suffix, std::cmp::min(forward_seq.len(), len_encoded))
}

pub fn common_suffix_length_fr(
    forward_seq: &[u8],
    encoded_to_reverse: &[u64],
    start_encoded: usize, // in bases (forward)
    end_encoded: usize,   // in bases (forward)
) -> usize {
    let len_encoded_seq = end_encoded - start_encoded;
    #[cfg(debug_assertions)]
    {
        // check the user passed start <= end (i.e. the start and end are correctly in forward)
        // TODO add those tests everywhere
        debug_assert!(start_encoded <= end_encoded);
        debug_assert!(encoded_to_reverse.len() >= (len_encoded_seq + 31) / 32);
    }
    let encoded_self_aligned_rigth = EncoderFromTheEndRightAligned::new(forward_seq);
    let revcomp_other_aligned_rigth =
        RevCompIterSRA::new(encoded_to_reverse, start_encoded, end_encoded);

    let size_suffix = common_suffix_length_for_iter(
        encoded_self_aligned_rigth,
        revcomp_other_aligned_rigth.into_iter(),
    );
    std::cmp::min(
        size_suffix,
        std::cmp::min(forward_seq.len(), len_encoded_seq),
    )
}

pub fn common_suffix_length_rf(
    seq_to_reverse: &[u8],
    encoded_forward_seq: &[u64],
    start: usize,
    end: usize,
) -> usize {
    let len_encoded_seq = end - start;
    debug_assert!(encoded_forward_seq.len() >= (len_encoded_seq + 31) / 32);
    let self_revcomp_encoded_aligned_right = RevCompEncoderSRA::new(seq_to_reverse);

    let encoded_other_aligned_rigth = FusedReverseIterator::new(encoded_forward_seq, start, end);

    let size_suffix = common_suffix_length_for_iter(
        self_revcomp_encoded_aligned_right,
        encoded_other_aligned_rigth,
    );
    std::cmp::min(
        size_suffix,
        std::cmp::min(seq_to_reverse.len(), len_encoded_seq),
    )
}

pub fn common_suffix_length_rr(
    seq_to_reverse: &[u8],
    encoded_to_reverse: &[u64],
    start: usize,
    end: usize,
) -> usize {
    let len_encoded_seq = end - start;
    debug_assert!(encoded_to_reverse.len() >= (len_encoded_seq + 31) / 32);
    let self_revcomp_encoded_aligned_right = RevCompEncoderSRA::new(seq_to_reverse);
    let revcomp_other_aligned_rigth = RevCompIterSRA::new(encoded_to_reverse, start, end);
    let size_suffix = common_suffix_length_for_iter(
        self_revcomp_encoded_aligned_right,
        revcomp_other_aligned_rigth,
    );
    std::cmp::min(
        size_suffix,
        std::cmp::min(seq_to_reverse.len(), len_encoded_seq),
    )
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::codec::{bmi::RevCompEncoder, Encoder};

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
        // TODO test with non 0 start
        {
            let a_string = String::from("ATCGGCGCATCG");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("ATCCGCGCATCG");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATGCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_ff(a, &b, 0, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_fr() {
        // TODO add more tests with start different than 0
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

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("ATCCGCGCATCG");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATGCGCGATCCGGCGCATCGGCAGTACGCATG");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_fr(a, &b, 0, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_rf() {
        // TODO test non 0
        {
            let a_string = String::from("CGATGCGCCGAT");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("CGATGCGCGGAT");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("TATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGCATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = Encoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rf(a, &b, 0, b_string.len()), 0);
        }
    }

    #[test]
    fn test_common_suffix_length_rr() {
        // TODO test non 0
        {
            let a_string = String::from("CGATGCGCCGAT");
            let b_string = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 12);
        }

        {
            //                                   ---v--------
            let a_string: String = String::from("CGATGCGCGGAT");
            let b_string: String = String::from("ATCGGCGCATCG");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 8);
        }
        {
            let a_string: String =
                String::from("ATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCAT");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 62);
        }

        {
            let a_string: String =
                String::from("TATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 63);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGAATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 0);
        }

        {
            let a_string: String =
                String::from("CATGCGTACTGCCGATGCGCCGGATCGCGCATGTCCTATAATCGTCGTCGTAGTCAGCGCTAT");
            let b_string: String =
                String::from("ATAGCGCTGACTACGACGACGATTATAGGACATTCGCGATCCGGCGCATCGGCAGTACGCATA");
            let a = a_string.as_bytes();
            let b = RevCompEncoder::new(b_string.as_bytes()).collect_vec();

            assert_eq!(common_suffix_length_rr(a, &b, 0, b_string.len()), 0);
        }
    }
}
