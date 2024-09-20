pub struct AlignLeftIterator<'a> {
    data: &'a [u64],
    start: usize,
    end: usize,
    current_pos: usize,
}

impl<'a> AlignLeftIterator<'a> {
    pub fn new(data: &'a [u64], start: usize, end: usize) -> Self {
        assert!(end <= data.len() * 32);
        AlignLeftIterator {
            data,
            start: start * 2,       // *2 to pass in bits
            end: end * 2,           // *2 to pass in bits
            current_pos: start * 2, // in in bits
        }
    }
}

impl<'a> Iterator for AlignLeftIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_pos >= self.end {
            return None; // we're done iterating.
        }

        let bit_len = self.end - self.current_pos; // number of bits that I can yield
        let current_idx = self.current_pos / 64;
        let end_idx = self.end / 64;
        let start_offset = self.current_pos % 64;

        let value = if start_offset == 0 {
            let end_offset = 64 - (self.end % 64);
            if current_idx == end_idx {
                // bits are within the same u64
                // OPTIMIZE use mask
                (self.data[current_idx] >> end_offset) << (end_offset)
            } else {
                // bits are across two u64 elements
                self.data[current_idx]
            }
        } else {
            let end_offset = 64 - (self.end % 64);
            if current_idx == end_idx {
                // bits are within the same u64
                (self.data[current_idx] >> end_offset) << (start_offset + end_offset)
            } else if end_offset == 64 {
                // TODO unit test
                // take all the data up to the end
                self.data[current_idx] << start_offset
            } else {
                // bits are across two u64 elements
                let first_part = self.data[current_idx] << start_offset;
                let second_part = if bit_len >= 64 {
                    self.data[current_idx + 1] >> (64 - start_offset)
                } else {
                    let second_part = self.data[current_idx + 1] >> end_offset;
                    let second_part = second_part << end_offset;

                    second_part >> (64 - start_offset)
                };
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
}
