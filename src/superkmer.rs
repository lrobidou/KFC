use crate::subsequence::NoBitPacked;
use crate::subsequence::Subsequence;
use crate::Minimizer;

use super::two_bits;
use std::iter::Map;
use std::iter::Rev;

const REVCOMP_TAB_CHAR: [char; 255] = {
    let mut tab = ['A'; 255];
    tab[b'A' as usize] = 'T';
    tab[b'T' as usize] = 'A';
    tab[b'C' as usize] = 'G';
    tab[b'G' as usize] = 'C';
    tab
};

pub fn reverse_complement(seq: &[u8]) -> String {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB_CHAR.get_unchecked(*base as usize) })
        .collect()
}

pub fn reverse_complement_ascii_to_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB_CHAR.get_unchecked(*base as usize) } as u8)
        .collect()
}

// TODO duplication: move to a module
const REVCOMP_TAB: [u8; 255] = {
    let mut tab = [0; 255];
    tab[b'A' as usize] = b'T';
    tab[b'T' as usize] = b'A';
    tab[b'C' as usize] = b'G';
    tab[b'G' as usize] = b'C';
    tab
};

// TODO "discuss" is there no copy here ? What is the cost of moving references ?
pub fn reverse_complement_no_copy(
    seq: impl DoubleEndedIterator<Item = u8>,
) -> Map<Rev<impl DoubleEndedIterator<Item = u8>>, fn(u8) -> u8> {
    seq.rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(base as usize) })
}

#[derive(Debug, Clone, Copy, PartialEq)]

pub struct Superkmer<'a> {
    minimizer: u64,
    start_mini: usize,
    end_mini: usize,
    pub superkmer: Subsequence<NoBitPacked<'a>>,
}

impl<'a> Superkmer<'a> {
    pub fn new(
        read: &'a [u8],
        start_mini: usize,
        end_mini: usize,
        start_sk: usize,
        end_sk: usize,
        same_orientation: bool,
    ) -> Self {
        let minizer_subsequence = &read[start_mini..end_mini];
        let minimizer = if same_orientation {
            two_bits::encode_minimizer(minizer_subsequence.iter().copied())
        } else {
            two_bits::encode_minimizer(reverse_complement_no_copy(
                minizer_subsequence.iter().copied(),
            ))
        };

        Self {
            minimizer,
            start_mini,
            end_mini,
            superkmer: Subsequence::<NoBitPacked>::new(read, start_sk, end_sk, same_orientation),
        }
    }

    pub fn hash_superkmer(&self) -> u64 {
        self.superkmer.hash()
    }

    pub fn get_minimizer(&self) -> Minimizer {
        self.minimizer
    }

    pub fn start_of_minimizer(&self) -> usize {
        self.start_mini
    }

    pub fn end_of_minimizer(&self) -> usize {
        self.end_mini
    }

    pub fn is_canonical_in_the_read(&self) -> bool {
        self.superkmer.same_orientation()
    }

    pub fn get_read(&self) -> &[u8] {
        self.superkmer.get_read()
    }

    #[cfg(debug_assertions)]
    pub fn minimizer_string(&self) -> String {
        use itertools::Itertools;
        if self.is_canonical_in_the_read() {
            String::from_utf8(
                self.get_read()[self.start_of_minimizer()..self.end_of_minimizer()]
                    .iter()
                    .copied()
                    .collect_vec(),
            )
            .unwrap()
        } else {
            reverse_complement(&self.get_read()[self.start_of_minimizer()..self.end_of_minimizer()])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let revcomp = reverse_complement(bytes);
        assert_eq!(revcomp, String::from("TAGCTCAGT"));
    }

    #[test]
    fn test_reverse_complement_ascii_to_ascii() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let revcomp = reverse_complement_ascii_to_ascii(bytes);
        assert_eq!(revcomp, String::from("TAGCTCAGT").as_bytes());
    }

    #[test]
    fn test_get_minimizer() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let sk = Superkmer::new(bytes, 1, 4, 1, 6, true);
        assert_eq!(sk.minimizer_string(), String::from("CTG"));
        let sk = Superkmer::new(bytes, 1, 4, 1, 6, false);
        assert_eq!(sk.minimizer_string(), String::from("CAG"));

        // TODO weird case: the minimizer is allowed to be outside of the superkmer
        let sk = Superkmer::new(bytes, 1, 4, 2, 6, false);
        assert_eq!(sk.minimizer_string(), String::from("CAG"));
    }
}
