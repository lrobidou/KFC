use super::superkmer::Superkmer;
use std::cmp::Ordering;
use std::iter::{Copied, Map, Rev};
use std::slice::Iter;

use crate::brrr_minimizers::MinimizerQueue;
// use crate::superkmer::SubsequenceMetadata;

// Get the reverse complement of a DNA sequence
pub fn reverse_complement<'a>(seq: &'a str) -> Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8> {
    seq.as_bytes().iter().rev().map(|base| match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => *base,
    })
}

pub fn same_orientation(seq: &str) -> Copied<Iter<'_, u8>> {
    seq.as_bytes().iter().copied()
}

pub fn is_canonical(seq: &str) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        if xc < yc {
            return true;
        } else if xc > yc {
            return false;
        }
    }
    // in case of palindrome, prefer saying the sequence is canonical
    true
}

pub enum SequenceOrientation<'a> {
    Same(Copied<Iter<'a, u8>>),
    Reverse(Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8>),
}

// a sequence in canonical form
// stores the sequence as it appears in the read and a boolean to read it backward if necessary
// TODO Clone is cheap, isn't it ?
#[derive(PartialEq, Eq, Clone, Debug, Copy)]
pub struct OrientedSequence<'a> {
    pub sequence: &'a str,
    pub is_same_orientation: bool,
}

impl<'a> OrientedSequence<'a> {
    pub fn new(sequence: &'a str, should_be_stored_in_the_same_orientation: bool) -> Self {
        Self {
            sequence,
            is_same_orientation: should_be_stored_in_the_same_orientation,
        }
    }

    pub fn is_same_orientation(&self) -> bool {
        self.is_same_orientation
    }

    pub fn get_oriented_sequence(&self) -> SequenceOrientation {
        if self.is_same_orientation {
            SequenceOrientation::Same(same_orientation(self.sequence))
        } else {
            SequenceOrientation::Reverse(reverse_complement(self.sequence))
        }
    }
}

fn compare(x: &mut impl Iterator<Item = u8>, y: &mut impl Iterator<Item = u8>) -> Ordering {
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        if xc < yc {
            return Ordering::Less;
        } else if xc > yc {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

impl<'a> PartialOrd for OrientedSequence<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for OrientedSequence<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        use SequenceOrientation as SO;
        match (self.get_oriented_sequence(), other.get_oriented_sequence()) {
            (SO::Same(mut x), SO::Same(mut y)) => compare(&mut x, &mut y),
            (SO::Same(mut x), SO::Reverse(mut y)) => compare(&mut x, &mut y),
            (SO::Reverse(mut x), SO::Same(mut y)) => compare(&mut x, &mut y),
            (SO::Reverse(mut x), SO::Reverse(mut y)) => compare(&mut x, &mut y),
        }
    }
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
struct MinimizerInfos<'a> {
    pub mmer: OrientedSequence<'a>,
    pub position: usize,
}

impl<'a> PartialOrd for MinimizerInfos<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> Ord for MinimizerInfos<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mmer.cmp(&other.mmer)
    }
}

// Function to compute superkmers
pub fn compute_superkmers_linear(sequence: &str, k: usize, m: usize) -> Vec<Superkmer> {
    let mut superkmers = Vec::with_capacity(sequence.len());
    if sequence.len() < k {
        return superkmers;
    }

    let w = k - m + 1;

    let mut mmers = Vec::with_capacity(sequence.len());
    for position in 0..sequence.len() - m + 1 {
        let mmer = OrientedSequence::new(
            &sequence[position..m + position],
            is_canonical(&sequence[position..m + position]),
        );
        mmers.push(MinimizerInfos { mmer, position });
    }

    let mut queue = MinimizerQueue::<_>::new(w);
    let mut minimizers: Vec<MinimizerInfos> = vec![];

    for (i, input) in mmers.iter().enumerate() {
        if i < w {
            queue.insert(input);
        } else {
            minimizers.push(*queue.get_min());

            queue.insert(input);
        }
    }
    minimizers.push(*queue.get_min());

    let mut current_minimizer = minimizers[0];
    let mut start_of_superkmer_as_read = 0;
    let mut end_of_superkmer_as_read = k;
    for (i, candidate_minimizer) in minimizers
        .iter()
        .enumerate()
        .take((sequence.len() - k) + 1)
        .skip(1)
    {
        // does the minimizer of this kmer have the same same position as the minimizer of the previous kmer?
        if current_minimizer.position == candidate_minimizer.position {
            end_of_superkmer_as_read += 1; // extend superkmer and keep other data untouched
        } else {
            superkmers.push(Superkmer::new(
                sequence,
                current_minimizer.position,
                current_minimizer.position + m,
                start_of_superkmer_as_read,
                end_of_superkmer_as_read,
                current_minimizer.mmer.is_same_orientation(),
            ));
            start_of_superkmer_as_read = i;
            end_of_superkmer_as_read = i + k;
            current_minimizer = *candidate_minimizer;
        }
    }

    // Add the last superkmer
    superkmers.push(Superkmer::new(
        sequence,
        current_minimizer.position,
        current_minimizer.position + m,
        start_of_superkmer_as_read,
        end_of_superkmer_as_read,
        current_minimizer.mmer.is_same_orientation(),
    ));

    superkmers
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTGCA").collect();
        assert_eq!(revcomp, b"TGCACTGCACAGT");
    }

    #[test]
    fn test_reverse_complement_n() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTNNGNCA").collect();
        assert_eq!(revcomp, b"TGNCNNACTGCACAGT");
    }

    #[test]
    fn test_get_canonical_kmer() {
        let kmer = "ACTGCGATGACGCAGATAGCAGATAGC";
        let canonical = OrientedSequence::new(kmer, is_canonical(kmer));

        assert!(canonical.is_same_orientation());
        match canonical.get_oriented_sequence() {
            SequenceOrientation::Same(seq_iter) => {
                assert_eq!(
                    seq_iter.collect::<Vec<u8>>(),
                    b"ACTGCGATGACGCAGATAGCAGATAGC"
                );
            }
            SequenceOrientation::Reverse(_) => {
                unreachable!();
            }
        }
    }

    #[test]
    fn test_get_canonical_kmer2() {
        let kmer = "GCTATCTGCTATCTGCGTCATCGCAGT";
        let canonical = OrientedSequence::new(kmer, is_canonical(kmer));
        assert!(!canonical.is_same_orientation());
        match canonical.get_oriented_sequence() {
            SequenceOrientation::Same(_) => {
                unreachable!();
            }
            SequenceOrientation::Reverse(seq_iter) => {
                assert_eq!(
                    seq_iter.collect::<Vec<u8>>(),
                    b"ACTGCGATGACGCAGATAGCAGATAGC"
                );
            }
        }
    }

    #[test]
    fn test_compute_superkmers() {
        let seq = "AGCAGCTAGCATTTTTGCAGT";
        assert_eq!(
            // ACTGCAAAAATGCTAGCTGCT
            compute_superkmers_linear(seq, 16, 5),
            vec![Superkmer::new(seq, 11, 11 + 5, 0, seq.len(), false)]
        );
    }
    // encoding and decoding
    #[test]
    fn test_compute_superkmers2() {
        // same minimizer further away (AAAAA)
        let seq = "AGCAGCTAGCATTTTTGCAGAAAAACC";
        assert_eq!(
            compute_superkmers_linear(seq, 16, 5),
            // AGCAGCTAGCATTTTTGCAGAAAA"
            //          CATTTTTGCAGAAAAACC
            vec![
                Superkmer::new(seq, 11, 11 + 5, 0, 24, false),
                Superkmer::new(seq, 20, 20 + 5, 9, 27, true),
            ]
        );
    }

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
        assert_eq!(compute_superkmers_linear("AGCAGCTAGCATTTT", 16, 5), vec![]);
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
}
