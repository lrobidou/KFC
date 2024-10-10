//          complete u64                      complete u64         incomplete u64
//  <------------------------------> <------------------------------> <------->
// |--------------------------------|--------------------------------|---------
//
// we put the incomlete u64 in a cache and fill it with a part of a complete u64
// we will return that
// but before we clean the cache and the rest of the complete u64 goes in the cache
pub struct FusedReverseIterator<'a> {
    // data: std::slice::Iter<'a, u64>,
    data: &'a [u64],
    start: usize, // in bits, in forward (so we will get closer to that at each call to `next``)
    current_pos: usize, // in bits
    nb_bits_in_cache: usize, // in [1:64]
    cache: u64,   // leftover bits from the previous iteration (aligned to the left)
    start_idx: usize, // index of the u64 in which start belongs
}

impl<'a> FusedReverseIterator<'a> {
    pub fn new(data: &'a [u64], start_in_forward: usize, end_in_forward: usize) -> Self {
        let start = start_in_forward * 2; //in bits
        let end = end_in_forward * 2; // in bits
        let end_idx = if end == 0 { 0 } else { (end - 1) / 64 };
        let nb_bits_in_cache = if end % 64 == 0 { 64 } else { end % 64 };

        // let data = &data[0..end_idx];
        // let mut data: std::slice::Iter<'_, u64> = data.iter();
        // let end_u64 = data.next().unwrap_or(&0);
        let end_u64 = data[end_idx];

        let cache = if nb_bits_in_cache == 64 {
            end_u64
        } else {
            end_u64 >> (64 - nb_bits_in_cache)
        };

        // position of the last base in the vector of u64
        let start_idx = if end_idx == 0 { 0 } else { (end_idx - 1) / 64 };

        FusedReverseIterator {
            data,
            cache,
            nb_bits_in_cache,
            start,
            current_pos: end,
            start_idx,
        }
    }
}

impl<'a> Iterator for FusedReverseIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_pos <= self.start {
            return None; // we're done iterating
        }

        debug_assert_eq!(
            if self.current_pos % 64 == 0 {
                64
            } else {
                self.current_pos % 64
            },
            self.nb_bits_in_cache
        );

        let bit_len = self.current_pos - self.start; // number of bits that I can yield

        let (value, cache, nb_advance) = if self.current_pos % 64 == 0 {
            let current_idx = (self.current_pos / 64) - 1;
            if self.start_idx < current_idx {
                // there is something before
                // safe as there is something before
                let previous_u64 = self.data[current_idx - 1];
                (self.cache, previous_u64, 64usize)
            } else {
                // we cannot assume there is anything before
                debug_assert!(self.start_idx == current_idx);
                // clear the low bits of the cache
                let val = (self.cache << (self.start % 64)) >> (self.start % 64);
                // cache is 0 as there is nothing after this iteration
                (val, 0, bit_len)
            }
        } else {
            // current pos is in the middle of a u64
            let current_idx = self.current_pos / 64;

            if self.start_idx < current_idx {
                // there is something before
                // safe as there is something before
                let previous_u64 = self.data[current_idx - 1];

                // let's complete the cache
                let prev = previous_u64;
                let down = prev << (self.current_pos % 64);
                let up = prev >> (64 - (self.current_pos % 64));
                let complete_u64 = down | self.cache;
                self.cache = up;
                if bit_len >= 64 {
                    // 64 bits are to be drawn, 64 of them will be returned
                    (complete_u64, up, 64)
                } else {
                    let nb_bits_to_throw = 64 - bit_len;
                    let val = (complete_u64 << nb_bits_to_throw) >> nb_bits_to_throw;
                    // returning a partial u64
                    // cache is 0 as there is nothing after this iteration

                    (val, 0, bit_len)
                }
            } else {
                // the start position and the current position are in the same index

                debug_assert!(self.start_idx == current_idx);
                let nb_bits_to_throw = 64 - bit_len;
                let val = (self.cache << nb_bits_to_throw) >> nb_bits_to_throw;
                (val, 0, bit_len)
            }
        };
        self.cache = cache;
        self.current_pos -= nb_advance;

        Some(value)
    }
}

pub fn revcomp_32_bases(encoded: u64) -> u64 {
    // reverse the bits
    // exchange odd and even bits
    // then flip odd bits to complement
    let reverse = encoded.reverse_bits();
    let odd: u64 = reverse & 0xAAAAAAAAAAAAAAAA;
    let even: u64 = reverse & 0x5555555555555555;
    let reverse = (odd >> 1) | (even << 1);
    reverse ^ 0xAAAAAAAAAAAAAAAA
}

pub fn revcomp_up_to_32_bases(encoded: u64, nb_bases: usize) -> u64 {
    let overful_revcomp = revcomp_32_bases(encoded);
    overful_revcomp << ((32 - nb_bases) * 2)
}

pub fn revcomp_up_to_32_bases_right_aligned(encoded: u64) -> u64 {
    revcomp_32_bases(encoded)
}

/// Takes a read encoded forward and iterates over the encoded reverse complement
pub struct RevCompIter<'a> {
    nb_bases_left: usize,
    reverse_iterator: FusedReverseIterator<'a>,
}

impl<'a> RevCompIter<'a> {
    pub fn new(data: &'a [u64], start_in_forward: usize, end_in_forward: usize) -> Self {
        let reverse_iterator = FusedReverseIterator::new(data, start_in_forward, end_in_forward);
        Self {
            reverse_iterator,
            nb_bases_left: end_in_forward - start_in_forward,
        }
    }
}

impl<'a> Iterator for RevCompIter<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        let next_element = self.reverse_iterator.next()?;
        debug_assert_ne!(self.nb_bases_left, 0);
        let r = if self.nb_bases_left >= 32 {
            Some(revcomp_32_bases(next_element))
        } else {
            // the next element is --------------------xxxxxxxxxx
            // where - is 0 and xxx are bases
            // let's correct that
            let left_align_element = next_element << ((32 - self.nb_bases_left) * 2);
            Some(revcomp_up_to_32_bases(
                left_align_element,
                self.nb_bases_left,
            ))
        };
        self.nb_bases_left = if self.nb_bases_left > 32 {
            self.nb_bases_left - 32
        } else {
            0
        };
        r
    }
}

/// Takes a read encoded forward and iterates over the encoded reverse complement, starting from the suffix.
/// The returned bases are aligned to the right.
pub struct RevCompIterSRA<'a> {
    data: &'a [u64],    // bases encoded
    end: usize,         // in bits, in forward
    current_pos: usize, // in bits, in forward
}

impl<'a> RevCompIterSRA<'a> {
    pub fn new(data: &'a [u64], start_in_forward: usize, end_in_forward: usize) -> Self {
        Self {
            data,
            end: end_in_forward * 2,           // pass in bits
            current_pos: start_in_forward * 2, // start from the start in forward, as it contains the suffix of the reverse complement
        }
    }
}

impl<'a> Iterator for RevCompIterSRA<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        // strategy: we group bases 64 by 64 for as long as there are 64 bases left
        // if not, we merge the rest
        // for each group, we pass in revcomp
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

                let merged_bases = self.data[current_idx] & mask;
                revcomp_up_to_32_bases_right_aligned(merged_bases)
            } else {
                // bits are in a single u64 already
                revcomp_32_bases(self.data[current_idx])
            }
        } else {
            // the start is not at the beginning of a u64
            if current_idx == end_idx {
                let data_to_remove_at_the_end = 64 - (bit_len + start_offset);
                debug_assert!(64 > data_to_remove_at_the_end);
                // bits are within the same u64
                let merged_bases = (self.data[current_idx] >> data_to_remove_at_the_end)
                    << (start_offset + data_to_remove_at_the_end);
                revcomp_up_to_32_bases_right_aligned(merged_bases)
            } else {
                // start offset != 0
                // spanning two u64
                // bits are across two u64 elements
                // (these two elements are valid)
                let first_part = self.data[current_idx] << start_offset;
                let second_part = self.data[current_idx + 1] >> (64 - start_offset);
                let merged_bases = first_part | second_part;
                // OPTIMIZE: we can cache data here
                // DEBUG check it works at the end of the read
                // (it may work for computing the suffix, but may add garbage at the end of the revcomp)
                revcomp_up_to_32_bases_right_aligned(merged_bases)
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

    use crate::codec::{Decoder, Encoder};

    use super::*;

    #[test]
    fn test_revcomp_up_to_32_bases() {
        let s = String::from("ATGACGACTGCTTACTCGATG");
        let s = s.as_bytes();

        let encoder = Encoder::new(s);
        let encoded = encoder.collect_vec()[0];
        let revcomp = revcomp_up_to_32_bases(encoded, s.len());
        assert_eq!(
            revcomp,
            // C A T C  G A G T  A A G C  A G T C  G T C A  T
            0b01001001_11001110_00001101_00111001_11100100_10000000_00000000_00000000
        );
    }

    #[test]
    fn test_fused_reverse_iterator_incomplete_bits() {
        let data = vec![
            0b1111000011110000111100001111000011110000111100001111000011110000u64,
            0b0000000000000000000000000000000000000000000000001111111111111111u64,
            0b1111111111111111000000000000000000000000000000000000000000000000u64,
        ];
        let data_rev = [
            0b0000000000000000000000000000000000001111111111111111111111111111u64,
            0b0000111100001111000011110000111100001111000011110000000000000000u64,
            0b0000000000000000000000000000000000000000000000000000111100001111u64,
        ];

        // suppose we have 70 pairs (140 bits) -> so the last u64 has only 12 bits of valid data.
        let num_pairs = 70;

        let iter = FusedReverseIterator::new(&data, 0, num_pairs).collect_vec();

        assert_eq!(iter, data_rev)
    }

    #[test]
    fn test_fused_reverse_iterator_incomplete_bits_start_middle() {
        let data = vec![
            0b1111000011110000111100001111000011110000111100001111000011110000u64,
            0b0000000000000000000000000000000000000000000000001111111111111111u64,
            0b1111111111111111000000000000000000000000000000000000000000000000u64,
        ];
        let data_rev = [
            0b0000000000000000000000000000000000001111111111111111111111111111u64,
            0b0000111100001111000011110000111100001111000011110000000000000000u64,
            0b0000000000000000000000000000000000000000000000000000000000000011u64,
        ];

        // suppose we start at 5 and end at 70 pairs (140 bits) -> so the last u64 has only 12 bits of valid data.
        let start = 5;
        let end = 70;

        let iter = FusedReverseIterator::new(&data, start, end).collect_vec();

        assert_eq!(iter, data_rev)
    }

    #[test]
    fn test_fused_reverse_iterator_one_incomplete_bits() {
        let data = vec![0b1111111100111111000000000000000000000000000000000000000000000000u64];
        let data_rev = [0b0000000000000000000000000000000000000000000000000000001111111100u64];
        let num_pairs = 5;

        let iter = FusedReverseIterator::new(&data, 0, num_pairs).collect_vec();

        assert_eq!(iter, data_rev)
    }

    #[test]
    fn test_fused_reverse_iterator_12_bases() {
        let read = String::from("ATCGGCGCATCG");
        let revcomp_read = String::from("CGATGCGCCGAT");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            assert_eq!(
                encoded,
                //  A T C G G C G C A T C G
                [0b0010011111011101001001110000000000000000000000000000000000000000]
            );
            let revcomp = RevCompIter::new(&encoded, 0, bytes.len()).collect_vec();
            assert_eq!(
                revcomp,
                //  C G A T G C G C C G A T
                [0b0111001011011101011100100000000000000000000000000000000000000000]
            );
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }

    #[test]
    fn test_fused_reverse_iterator_32_bases() {
        let read = String::from("GTTCCCGAGTACTTCACTTATACTAGACATGG");
        let revcomp_read = String::from("CCATGTCTAGTATAAGTGAAGTACTCGGGAAC");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            assert_eq!(
                encoded,
                //  G T T C C C G A G T A C T T C A C T T A T A C T A G A C A T G G
                [0b1110100101011100111000011010010001101000100001100011000100101111]
            );
            let revcomp = RevCompIter::new(&encoded, 0, bytes.len()).collect_vec();
            assert_eq!(
                revcomp,
                //  C C A T G T C T A G T A T A A G T G A A G T A C T C G G G A A C
                [0b0101001011100110001110001000001110110000111000011001111111000001]
            );
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }

    #[test]
    fn test_fused_reverse_iterator_64_bases() {
        let read = String::from("GGACTTGTTTGCCTCTTTGTAATACGTAGGAGCTTGATGATAAGTCTAAGTAATCACATGATCG");
        let revcomp_read =
            String::from("CGATCATGTGATTACTTAGACTTATCATCAAGCTCCTACGTATTACAAAGAGGCAAACAAGTCC");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            let revcomp = RevCompIter::new(&encoded, 0, bytes.len()).collect_vec();
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }

    #[test]
    fn test_fused_reverse_iterator_100_bases() {
        let read = String::from("AAGTTAGGGGCTGTGTGCAAAAGTTAGGGGCGGTATTAGCTTGAATAAGTGCTCTACCTCGACCTTATCCGACCATGAATGGCCACGGTCAGGCGAGTGG");
        let revcomp_read =
            String::from("CCACTCGCCTGACCGTGGCCATTCATGGTCGGATAAGGTCGAGGTAGAGCACTTATTCAAGCTAATACCGCCCCTAACTTTTGCACACAGCCCCTAACTT");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            let revcomp = RevCompIter::new(&encoded, 0, bytes.len()).collect_vec();
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }

    //     input:
    // data  = 0101001100011100000100110001110101111011111000000101001111000100
    // data  = 1001000101010001001010000010001011011111011011000010001101110101
    // data  = 0111101011100001011100100001110100111100001111101101110000101101
    // data  = 0101000000000000000000000000000000000000000000000000000000000000
    // start_in_forward = 15
    // end_in_forward = 98
    // input:
    // data  = 0001110110111101110111101111011101001110010000111010111001110000
    // data  = 0001001100000011101100100111100100000001001001101111001011010100
    // data  = 1010010011000001101010111101110100001000000000000000000000000000
    // thread '<unnamed>' panicked at src/subsequence.rs:149:13:
    // assertion `left == right` failed
    //   left: 74
    //  right: 0

    #[test]
    fn test_recvomp_iter_sra() {
        let start = 15;
        let end = 98;
        let data = [
            /////////////////////////////// start of the first base //////////
            /////////////////////////////// v ////////////////////////////////
            0b0101001100011100000100110001110101111011111000000101001111000100,
            0b1001000101010001001010000010001011011111011011000010001101110101,
            0b0111101011100001011100100001110100111100001111101101110000101101,
            //   end of the last base ////////////////////////////////////////
            // / v ///////////////////////////////////////////////////////////
            0b0101000000000000000000000000000000000000000000000000000000000000,
        ];

        // suffix right aligned:
        //  0b0101111011111000000101001111000100100100010101000100101000001000
        //  = CCGTGGTAACCAGGACATCACCCACATTAATA
        // rc=TATTAATGTGGGTGATGTCCTGGTTACCACGG
        //  0b1011011111011011000010001101110101011110101110000101110010000111
        //  = TGCGGCTGAATAGCGCCCGTTGTACCGATACG
        // rc=CGTATCGGTACAACGGGCGCTATTCAGCCGCA
        //  0b0100111100001111101101110000101101010100000000000000000000000000
        //  = CAGGAAGGTGCGAATGCCCAAAAAAAAAAAAA
        // rc=TTTTTTTTTTTTTGGGCATTCGCACCTTCCTG

        let expected = vec![
            // T A T T A A T G T G G G T G A T G T C C T G G T T A C C A C G G
            0b1000101000001011101111111011001011100101101111101000010100011111,
            // C G T A T C G G T A C A A C G G G C G C T A T T C A G C C G C A
            0b0111100010011111100001000001111111011101100010100100110101110100,
            // T T T T T T T T T T T T T G G G C A T T C G C A C C T T C C T G
            0b1010101010101010101010101011111101001010011101000101101001011011,
        ];

        let result = RevCompIterSRA::new(&data, start, end).collect_vec();
        assert_eq!(expected, result);
    }
}
