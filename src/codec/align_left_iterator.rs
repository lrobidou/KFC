pub struct AlignLeftIterator<'a> {
    data: &'a [u64],
    end: usize,
    current_pos: usize,
}

impl<'a> AlignLeftIterator<'a> {
    pub fn new(data: &'a [u64], start: usize, end: usize) -> Self {
        assert!(end <= data.len() * 32);
        AlignLeftIterator {
            data,
            end: end * 2,           // *2 to pass in bits
            current_pos: start * 2, // *2 to pass in bits
        }
    }
}

impl Iterator for AlignLeftIterator<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_pos >= self.end {
            return None; // we're done iterating
        }

        let bit_len = self.end - self.current_pos; // number of bits that I can yield
        let current_idx = self.current_pos / 64;

        // position of the last base in the vector of u64
        let end_idx = if self.end == 0 {
            0
        } else {
            (self.end - 1) / 64
        };

        // number of bits I should ignore at the start of a u64
        let start_offset = self.current_pos % 64;

        let value = if start_offset == 0 {
            // are we at the last iteration, and is that last iteration not full ?
            if current_idx == end_idx && (self.end % 64) != 0 {
                let end_offset = 64 - (self.end % 64);
                // bits are within the same u64

                // equivalent to (self.data[current_idx] >> end_offset) << (end_offset)
                let mask = !((1 << end_offset) - 1);
                self.data[current_idx] & mask
            } else {
                // bits are in a single u64 already
                self.data[current_idx]
            }
        } else {
            // the start is not at the beginning of a u64
            if current_idx == end_idx {
                let data_to_remove_at_the_end = 64 - (bit_len + start_offset);
                debug_assert!(64 > data_to_remove_at_the_end);
                // bits are within the same u64
                (self.data[current_idx] >> data_to_remove_at_the_end)
                    << (start_offset + data_to_remove_at_the_end)
            } else {
                // start offset != 0
                // spanning two u64
                // bits are across two u64 elements
                // (these two elements are valid)
                let first_part = self.data[current_idx] << start_offset;
                let second_part = self.data[current_idx + 1] >> (64 - start_offset);
                first_part | second_part
            }
        };

        // Advance the current position
        self.current_pos += 64;

        Some(value)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::codec::Encoder;

    use super::*;

    #[test]
    fn test_align_left_iterator_end_same_u64() {
        let v: Vec<u64> = vec![0b1100010001010001100000001111001100101011000101010000011000000111];
        let extracted: Vec<u64> =
            vec![0b0001000101000110000000111100000000000000000000000000000000000000];
        let aligned = AlignLeftIterator::new(&v, 1, 14).collect_vec();
        assert_eq!(aligned, extracted);
    }

    #[test]
    fn test_align_left_iterator_end_different_u64() {
        let v: Vec<u64> = vec![
            0b1100010001010001100000001111001100101011000101010000011000000111,
            0b0001110001010001100000001111001101001101000101010000011000000111,
        ];
        let extracted: Vec<u64> = vec![
            0b0100011000000011110011001010110001010100000110000001110001110001,
            0b0100011000000011110011010011010001010100000110000000000000000000,
        ];
        let aligned = AlignLeftIterator::new(&v, 5, 60).collect_vec();
        assert_eq!(aligned, extracted);
    }

    #[test]
    fn test_all_g() {
        let read = String::from("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
        let encoded = Encoder::new(read.as_bytes()).collect_vec();
        let start = 34;
        let end = 128;
        let subpart = &read[start..end];
        let encoded_subpart = Encoder::new(subpart.as_bytes()).collect_vec();

        let aligned = AlignLeftIterator::new(&encoded, start, end).collect_vec();
        assert_eq!(aligned, encoded_subpart);
    }

    #[test]
    fn test_ko() {
        // let read = "AAATAATGGTGGCTACACGCATGAGAAATGAGAAGAAAAGAGGGTGATATATGTTTGAGTTTATTGGTTTTTTATTTTTGCTTTTGGTGTGTTATGTCTTTATACAAATGGTTGCTGTCAGGGTGTTTCCTGAATATGGAGTACGGAAAGAAAAAAAAAA";
        let read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
        let start = 113;
        let end = 150;
        let subseq = &read[start..end];

        let encoded = Encoder::new(read.as_bytes()).collect_vec();
        let encoded_subseq = Encoder::new(subseq.as_bytes()).collect_vec();
        let aligned = AlignLeftIterator::new(&encoded, start, end).collect_vec();
        assert_eq!(encoded_subseq, aligned);
        // [147474613200011443, 865605587707809442, 13738946082082057867, 11135221266529976042, 6559700535848468480], 113, 150
        // decoded = GCTGTCAGGGTGTTTCCTGAATATGGAGTAAAGAAAG
    }
}
