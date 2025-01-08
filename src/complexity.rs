//! Computes the complexity of k-mers in streaming

use crate::two_bits::char_to_2bit;

const W: u32 = 5; // length of the words used to compute the complexity of k-mers
const NB_W_MERS: usize = 4_usize.pow(W); // number of possible w-mers
const ENCODING_MASK: u16 = (1u16 << (2 * W)) - 1; // mask that selects the 2*w lowest bits

/// Add an ASCII base to the encoding aof a w-mer
fn add_base_to_encoding(encoding: u16, base: u8) -> u16 {
    let mut encoding = encoding;
    encoding <<= 2; // make place for the new base
    encoding &= ENCODING_MASK; // keep only relevant bits
    let encoded_base = char_to_2bit(base); //encode new base
    encoding += u16::from(encoded_base);
    encoding
}

/// Computes the complexity of k-mers in streaming
/// Complexity is defined as the number of distints [`W`]-mers in a `k`-mer.
pub struct ComplexityComputation<'a> {
    k: usize,                        // size of k-mers
    current_complexity: u16,         // complexity of the current k-mer
    superkmer: &'a [u8],             // superkmer from which the k-mers will be computed
    current_position: usize,         // position of the current k-mer in the superkmer
    current_encoding: u16,           // encoding of the last w-mer of the current k-mer
    old_encoding: u16,               // encoding of the first w-mer of the window
    current_counts: [u8; NB_W_MERS], // counts of w-mers in the window
}

impl<'a> ComplexityComputation<'a> {
    pub fn new(superkmer: &'a [u8], k: usize) -> Self {
        debug_assert!(superkmer.len() >= k);
        let mut current_encoding: u16 = 0;
        let mut current_complexity = 0;
        let mut current_counts = [0; 4_usize.pow(W)];

        // fill the encoding of the first w bases
        for base in superkmer.iter().take(W as usize) {
            current_encoding <<= 2;
            let encoded_base = char_to_2bit(*base);
            current_encoding += u16::from(encoded_base);
            debug_assert!(current_encoding == (current_encoding & ENCODING_MASK));
        }

        // this is the first w-mer
        let old_encoding = current_encoding;

        // we have our first w-mer, let's increase its count and the complexity
        current_counts[usize::from(current_encoding)] += 1;
        current_complexity += 1;

        // starting from here, any new base will for a valid w-mer
        // let's fill these w-mers up to one base before the first k-mer
        for base in superkmer.iter().take(k - 1).skip(W as usize) {
            current_encoding <<= 2; // make place for the new base
            current_encoding &= ENCODING_MASK; // keep only relevant bits
            let encoded_base = char_to_2bit(*base); //encode new base
            current_encoding += u16::from(encoded_base);
            if current_counts[usize::from(current_encoding)] == 0 {
                current_complexity += 1;
            }
            current_counts[usize::from(current_encoding)] += 1;
        }

        Self {
            k,
            superkmer,
            current_complexity,
            current_position: 0,
            old_encoding,
            current_encoding,
            current_counts,
        }
    }

    #[cfg(test)]
    fn get_counts(&self) -> [u8; NB_W_MERS] {
        self.current_counts
    }

    fn get_complexity_slow(&self) -> u16 {
        let mut complexity = 0;
        for count in self.current_counts {
            if count > 0 {
                complexity += 1;
            }
        }
        complexity
    }

    fn shift_old_window(&mut self, base: u8) {
        debug_assert!(self.current_counts[usize::from(self.old_encoding)] > 0);
        self.current_counts[usize::from(self.old_encoding)] -= 1;
        if self.current_counts[usize::from(self.old_encoding)] == 0 {
            self.current_complexity -= 1;
        }
        self.old_encoding = add_base_to_encoding(self.old_encoding, base);
    }

    fn shift_current_window(&mut self, base: u8) {
        self.current_encoding = add_base_to_encoding(self.current_encoding, base);
        if self.current_counts[usize::from(self.current_encoding)] == 0 {
            self.current_complexity += 1;
        }
        self.current_counts[usize::from(self.current_encoding)] += 1;
    }
}

impl Iterator for ComplexityComputation<'_> {
    type Item = u16;

    fn next(&mut self) -> Option<u16> {
        if self.current_position > self.superkmer.len() - self.k {
            return None;
        }

        if self.current_position > 0 {
            // we need to remove the first w-mer of the window
            // and update the encoding of the first w-mer of the window
            self.shift_old_window(self.superkmer[self.current_position + W as usize - 1]);
        }

        if self.current_position + self.k - 1 < self.superkmer.len() {
            self.shift_current_window(self.superkmer[self.current_position + self.k - 1]);
        }

        // the current k-mer's compliexity is known
        self.current_position += 1;

        debug_assert_eq!(self.get_complexity_slow(), self.current_complexity);

        Some(self.current_complexity)
    }
}

pub fn is_complexity_above_threshold(
    superkmer: &[u8],
    k: usize,
    complexity_threshold: u16,
) -> bool {
    if let Some(min_complexity) = ComplexityComputation::new(superkmer, k).max() {
        if min_complexity > complexity_threshold {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_complexity_one_kmer() {
        let superkmer = "AAAAAAAAAA".as_bytes();
        let k = 10;
        let mut counts = [0u8; NB_W_MERS];
        counts[0] = 5;
        let mut complexity_iterator = ComplexityComputation::new(superkmer, k);
        // start state
        assert_eq!(complexity_iterator.get_counts(), counts);
        // first iteration
        assert_eq!(complexity_iterator.next().unwrap(), 1);
        counts[0] = 6;
        assert_eq!(complexity_iterator.get_counts(), counts);
        //second iteration
        assert_eq!(complexity_iterator.next(), None);
    }

    #[test]
    fn test_complexity_two_kmer() {
        let superkmer = "AAAAAAAAAAA".as_bytes();
        let k = 10;

        let expected_complexity = vec![1, 1];

        let complexity_iterator = ComplexityComputation::new(superkmer, k);
        assert_eq!(complexity_iterator.collect_vec(), expected_complexity);
    }

    #[test]
    fn test_complexity_two_distinct_kmer() {
        let superkmer = "AAAAAAAAAAC".as_bytes();
        let k = 10;

        let expected_complexity = vec![1, 2];

        let complexity_iterator = ComplexityComputation::new(superkmer, k);
        assert_eq!(complexity_iterator.collect_vec(), expected_complexity);
    }

    #[test]
    fn test_complexity_5_distinct_wmers_0() {
        let superkmer = "ACGCTACGCTACGCTACGCT".as_bytes();
        let k = 20;

        let expected_complexity = vec![5; 1];

        let complexity_iterator = ComplexityComputation::new(superkmer, k);
        assert_eq!(complexity_iterator.collect_vec(), expected_complexity);
    }

    #[test]
    fn test_complexity_5_distinct_wmers_two_kmers() {
        let superkmer = "ACGCTACGCTACGCTACGCTA".as_bytes();
        let k = 20;

        let expected_complexity = vec![5; 2];

        let complexity_iterator = ComplexityComputation::new(superkmer, k);
        assert_eq!(complexity_iterator.collect_vec(), expected_complexity);
    }

    #[test]
    fn test_complexity_5_distinct_wmers_two_kmers_increase_complexity() {
        let superkmer = "ACGCTACGCTACGCTACGCTATGGCACGATTCGATT".as_bytes();
        //                     ^--- ACGCT will repeat
        //                           until here ---^
        //                                          from here, the complexity increases
        //                                          until it reaches the max (16)
        //                                               ^--- CGATT is reapeated from here
        //                                                    after one repetition, the complexity decreases
        let k = 20;

        let expected_complexity = vec![
            // repetition of ACGCT
            5, 5, // complexity increases
            6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, // complexity is the max possible
            16, 16, 16, 15, // complexity deacrese (repetition of CGATT)
        ];

        let complexity_iterator = ComplexityComputation::new(superkmer, k);
        assert_eq!(complexity_iterator.collect_vec(), expected_complexity);
    }

    #[test]
    fn random_kmers_are_not_low_complexity() {
        let superkmer = "TTCTACTGCTTAGTCATGGACGCCATTCGAGCCTATTGGGTCACTGTGCGCTGGGTTGAGTAGGCGCATGCGAAAACCTAAATGCCCTAACCATATTCCATTTCCGTATTTCCCCAACACTGGGTGGGCGTTATCTCGGGCCTATTTATTGTGAAGCGGGCTCGTCTTCTAACCTAATAAGGGTAGTGTTTGGTTCGATA".as_bytes();

        let complexity_20 = ComplexityComputation::new(superkmer, 20).min().unwrap();
        let complexity_30 = ComplexityComputation::new(superkmer, 30).min().unwrap();
        let complexity_50 = ComplexityComputation::new(superkmer, 50).min().unwrap();
        let complexity_100 = ComplexityComputation::new(superkmer, 100).min().unwrap();
        let complexity_200 = ComplexityComputation::new(superkmer, 200).min().unwrap();

        // not low complexity
        assert!(complexity_20 > 8);
        // k increases => complexity increases
        assert!(complexity_30 > complexity_20);
        assert!(complexity_50 > complexity_30);
        assert!(complexity_100 > complexity_50);
        assert!(complexity_200 > complexity_100);
    }
}
