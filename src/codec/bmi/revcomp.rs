//          complete u64                      complete u64         incomplete u64
//  <------------------------------> <------------------------------> <------->
// |--------------------------------|--------------------------------|---------
//
// we put the incomlete u64 in a cache and fill it with a part of a complete u64
// we will return that
// but before we clean the cache and the rest of the complete u64 goes in the cache
pub struct FusedReverseIterator<'a> {
    data: &'a [u64],
    current_bit: usize, // The bit index we're currently processing
    cache: u64,         // leftover bits from the previous iteration
    nb_bits_in_cache: usize,
}

impl<'a> FusedReverseIterator<'a> {
    pub fn new(data: &'a [u64], nb_bases: usize) -> Self {
        let nb_bits_in_cache = (nb_bases % 32) * 2;
        let cache = if nb_bits_in_cache == 0 {
            0
        } else {
            data[data.len() - 1] >> (64 - nb_bits_in_cache)
        };
        let current_bit = nb_bases * 2; // start from the last bit
        FusedReverseIterator {
            data,
            current_bit,
            cache,
            nb_bits_in_cache,
        }
    }
}

impl<'a> Iterator for FusedReverseIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_bit == 0 {
            return None; // no more bits to process
        }

        debug_assert_eq!(self.current_bit % 64, self.nb_bits_in_cache);

        if self.current_bit < 64 {
            self.current_bit = 0;
            return Some(self.cache);
        }
        // we want to gather 64 bits at a time
        let idx = (self.current_bit - self.nb_bits_in_cache) / 64 - 1;
        let val = self.data[idx];

        let val_low = val << self.nb_bits_in_cache; // make some space to merge in the cache
        let val_up = if self.nb_bits_in_cache == 0 {
            0 //no bits in cache => the cache is empty
        } else {
            val >> (64 - self.nb_bits_in_cache) // the remaining bits are the future cache
        };
        debug_assert_eq!(val_low & self.cache, 0); // no overlapp between the two parts I'm merging
        let r = val_low | self.cache;

        self.cache = val_up;
        self.current_bit -= 64;

        Some(r)
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

pub unsafe fn revcomp_up_to_32_bases(encoded: u64, nb_bases: usize) -> u64 {
    let overful_revcomp = revcomp_32_bases(encoded);
    overful_revcomp << ((32 - nb_bases) * 2)
}

/// Takes a read encoded forward and iterates over the encoded reverse complement
pub struct RevCompIter<'a> {
    nb_bases_left: usize,
    reverse_iterator: FusedReverseIterator<'a>,
}

impl<'a> RevCompIter<'a> {
    pub fn new(data: &'a [u64], nb_bases: usize) -> Self {
        let reverse_iterator = FusedReverseIterator::new(data, nb_bases);
        Self {
            reverse_iterator,
            nb_bases_left: nb_bases,
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
            Some(unsafe { revcomp_up_to_32_bases(left_align_element, self.nb_bases_left) })
        };
        self.nb_bases_left = if self.nb_bases_left > 32 {
            self.nb_bases_left - 32
        } else {
            0
        };
        r
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
        let revcomp = unsafe { revcomp_up_to_32_bases(encoded, s.len()) };
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

        let iter = FusedReverseIterator::new(&data, num_pairs).collect_vec();

        assert_eq!(iter, data_rev)
    }

    #[test]
    fn test_fused_reverse_iterator_one_incomplete_bits() {
        let data = vec![0b1111111100111111000000000000000000000000000000000000000000000000u64];
        let data_rev = [0b0000000000000000000000000000000000000000000000000000001111111100u64];
        let num_pairs = 5;

        let iter = FusedReverseIterator::new(&data, num_pairs).collect_vec();

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
            let revcomp = RevCompIter::new(&encoded, bytes.len()).collect_vec();
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
            let revcomp = RevCompIter::new(&encoded, bytes.len()).collect_vec();
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
            let revcomp = RevCompIter::new(&encoded, bytes.len()).collect_vec();
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
            let revcomp = RevCompIter::new(&encoded, bytes.len()).collect_vec();
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }
}
