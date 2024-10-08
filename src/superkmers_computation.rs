use super::superkmer::Superkmer;
use crate::minimizer_iter::{simd_minimizer_iter, MinimizerItem};
use crate::superkmer::REVCOMP_TAB;
use std::cmp::Ordering;
use std::iter::{Copied, Map, Rev};
use std::slice::Iter;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

// TODO should we do something else instead of using this complex type?
#[allow(clippy::type_complexity)]
fn reverse_complement<'a>(seq: &'a [u8]) -> Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8> {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(*base as usize) })
}

pub fn same_orientation(seq: &[u8]) -> Copied<Iter<'_, u8>> {
    seq.iter().copied()
}

pub fn is_canonical(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };

        match xc.cmp(&yc) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => {}
        }
    }
    // in case of palindrome, prefer saying the sequence is canonical
    true
}

pub struct SuperkmerIterator<'a, I: Iterator<Item = MinimizerItem>> {
    sequence: &'a [u8],
    minimizer_iter: std::iter::Peekable<I>,
    k: usize,
    m: usize,
    previous_minimizer: Option<MinimizerItem>,
}

impl<'a, I: Iterator<Item = MinimizerItem>> SuperkmerIterator<'a, I> {
    fn new(sequence: &'a [u8], minimizer_iter: I, k: usize, m: usize) -> Self {
        let mut minimizer_iter = minimizer_iter.peekable();
        let previous_minimizer: Option<MinimizerItem> = minimizer_iter.next();
        Self {
            sequence,
            minimizer_iter,
            k,
            m,
            previous_minimizer,
        }
    }
}

impl<'a, I: Iterator<Item = MinimizerItem>> Iterator for SuperkmerIterator<'a, I> {
    type Item = Superkmer<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.previous_minimizer {
            Some(ref previous_minimizer) => match self.minimizer_iter.next() {
                Some(current_minimizer) => {
                    let mid = previous_minimizer.position + self.m / 2;
                    let same_orientation = self.sequence[mid] & 0b100 == 0;
                    let sk = Superkmer::new(
                        self.sequence,
                        previous_minimizer.position,
                        previous_minimizer.position + self.m,
                        previous_minimizer.window_start,
                        current_minimizer.window_start + self.k - 1,
                        same_orientation,
                    );
                    self.previous_minimizer = Some(current_minimizer);
                    Some(sk)
                }
                None => {
                    let mid = previous_minimizer.position + self.m / 2;
                    let same_orientation = self.sequence[mid] & 0b100 == 0;
                    let sk = Superkmer::new(
                        self.sequence,
                        previous_minimizer.position,
                        previous_minimizer.position + self.m,
                        previous_minimizer.window_start,
                        self.sequence.len(),
                        same_orientation,
                    );
                    self.previous_minimizer = None;
                    Some(sk)
                }
            },
            None => None,
        }
    }
}

pub fn compute_superkmers_linear_streaming(
    sequence: &[u8],
    k: usize,
    m: usize,
) -> Option<SuperkmerIterator<impl Iterator<Item = MinimizerItem> + '_>> {
    if sequence.len() < k {
        None
    } else {
        let minimizer_iter = simd_minimizer_iter(sequence, m, k - m + 1);
        let superkmer_iter = SuperkmerIterator::new(sequence, minimizer_iter, k, m);
        Some(superkmer_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTGCA".as_bytes()).collect();
        assert_eq!(revcomp, b"TGCACTGCACAGT");
    }

    #[test]
    fn test_reverse_complement_n() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTNNGNCA".as_bytes()).collect();
        // assert_eq!(revcomp, b"TG\0C\0\0ACTGCACAGT");
        assert_eq!(revcomp, b"TGTCTTACTGCACAGT");
    }

    // #[test]
    // fn test_get_canonical_kmer() {
    //     let kmer = "ACTGCGATGACGCAGATAGCAGATAGC".as_bytes();
    //     let canonical = OrientedSequence::new(kmer, is_canonical(kmer));

    //     assert!(canonical.is_same_orientation());
    //     match canonical.get_oriented_sequence() {
    //         SequenceOrientation::Same(seq_iter) => {
    //             assert_eq!(
    //                 seq_iter.collect::<Vec<u8>>(),
    //                 b"ACTGCGATGACGCAGATAGCAGATAGC"
    //             );
    //         }
    //         SequenceOrientation::Reverse(_) => {
    //             unreachable!();
    //         }
    //     }
    // }

    // #[test]
    // fn test_get_canonical_kmer2() {
    //     let kmer = "GCTATCTGCTATCTGCGTCATCGCAGT".as_bytes();
    //     let canonical = OrientedSequence::new(kmer, is_canonical(kmer));
    //     assert!(!canonical.is_same_orientation());
    //     match canonical.get_oriented_sequence() {
    //         SequenceOrientation::Same(_) => {
    //             unreachable!();
    //         }
    //         SequenceOrientation::Reverse(seq_iter) => {
    //             assert_eq!(
    //                 seq_iter.collect::<Vec<u8>>(),
    //                 b"ACTGCGATGACGCAGATAGCAGATAGC"
    //             );
    //         }
    //     }
    // }
    // encoding[b'A' as usize] = 0b00;
    // encoding[b'C' as usize] = 0b01;
    // encoding[b'G' as usize] = 0b10;
    // encoding[b'T' as usize] = 0b11;
    // 11 01 10 00 11
    // TCGAT
    // ATCGA
    // #[test]
    // fn test_compute_superkmers() {
    //     let seq = "AGCAGCTAGCATTTTTGCAGT".as_bytes();
    //     let superkmers: Vec<Superkmer> = compute_superkmers_linear_streaming(seq, 17, 5)
    //         .unwrap()
    //         .collect();
    //     assert_eq!(
    //         // ACTGCAAAAATGCTAGCTGCT
    //         superkmers,
    //         vec![Superkmer::new(seq, 11, 11 + 5, 0, seq.len(), false)]
    //     );
    // }
    // encoding and decoding
    // TODO correct test
    // #[test]
    // fn test_compute_superkmers2() {
    //     // same minimizer further away (AAAAA)
    //     let seq = "AGCAGCTAGCATTTTTGCAGAAAAACC".as_bytes();
    //     let superkmers: Vec<Superkmer> = compute_superkmers_linear_streaming(seq, 17, 5)
    //         .unwrap()
    //         .collect();
    //     assert_eq!(
    //         superkmers,
    //         // AGCAGCTAGCATTTTTGCAGAAAA
    //         //         GCATTTTTGCAGAAAAACC
    //         // AGCAGCTAGCATTTTTGCAGAAAAAC
    //         //           ATTTTTGCAGAAAAACC
    //         vec![
    //             Superkmer::new(seq, 11, 11 + 5, 0, 24, false),
    //             Superkmer::new(seq, 20, 20 + 5, 10, 27, true),
    //         ]
    //     );
    // }

    // #[test]
    // fn test_compute_superkmers3() {
    //     let superkmers = compute_superkmers_linear("AAGACGCGCCAGCGTCGCATCAGGCGTTGAATGCCGGATGCGCTTCCTGATAAGACGCGCCAGCGTCGCATCAGGCGTTGAATGCCGGATGCGCTT", 40, 31);
    //     for x in &superkmers {
    //         println!("{:?}", x);
    //     }
    //     assert_eq!(superkmers, vec![]);
    // }

    #[test]
    fn test_compute_superkmers_sequence_too_short() {
        // there are no superkmers for sequence < k
        if let Some(_x) = compute_superkmers_linear_streaming("AGCAGCTAGCATTTT".as_bytes(), 16, 5) {
            panic!()
        }
    }

    // #[test]
    // fn test_print_() {
    //     let super_kmer_infos = SuperKmerInfos {
    //         superkmer: String::from("CATTTTTGCAGAAAAACC"),
    //         minimizer: String::from("AAAAA"),
    //         was_read_canonical: true,
    //         start_of_minimizer_as_read: 20,
    //         start_of_superkmer_as_read: 9,
    //     };
    //     let s = "SuperKmerInfos { superkmer: \"CATTTTTGCAGAAAAACC\", minimizer: \"AAAAA\", was_read_canonical: true, start_of_minimizer_as_read: 20, start_of_superkmer_as_read: 9 }";
    //     assert_eq!(format!("{:?}", super_kmer_infos), s);
    // }
    // #[test]
    // fn test_minimizer_iter() {
    //     let seq = b"TGATTGC";
    //     let minimizer_size = 3;
    //     let width = 4;
    //     let hasher = BuildNoHashHasher::<u64>::default();
    //     let mut min_iter = MinimizerBuilder::new()
    //         .minimizer_size(minimizer_size)
    //         .width(width)
    //         .hasher(hasher)
    //         .iter(seq);

    //     assert_eq!(min_iter.next(), Some((0b001111, 2))); // ATT
    //                                                       // assert_eq!(min_iter.next(), Some((0b010001, 6))); // CAC
    //                                                       // assert_eq!(min_iter.next(), Some((0b000100, 7))); // ACA
    //                                                       // assert_eq!(min_iter.next(), Some((0b000011, 9))); // AAT
    //     assert_eq!(min_iter.next(), None);
    // }
}
