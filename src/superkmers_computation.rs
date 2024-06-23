use super::superkmer::Superkmer;
use std::cmp::Ordering;
use std::iter::{Copied, Map, Rev};
use std::slice::Iter;

use crate::brrr_minimizers::MinimizerQueue;
// use crate::superkmer::SubsequenceMetadata;

// Get the reverse complement of a DNA sequence
pub fn reverse_complement<'a>(seq: &'a [u8]) -> Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8> {
    seq.iter().rev().map(|base| match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => *base,
    })
}

pub fn same_orientation(seq: &[u8]) -> Copied<Iter<'_, u8>> {
    seq.iter().copied()
}

pub fn is_canonical(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        match xc.cmp(&yc) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => {}
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
    pub sequence: &'a [u8],
    pub is_same_orientation: bool,
}

impl<'a> OrientedSequence<'a> {
    pub fn new(sequence: &'a [u8], should_be_stored_in_the_same_orientation: bool) -> Self {
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
        match xc.cmp(&yc) {
            Ordering::Less => return Ordering::Less,
            Ordering::Greater => return Ordering::Greater,
            Ordering::Equal => {}
        }
    }
    Ordering::Equal
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

impl<'a> PartialOrd for OrientedSequence<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct MmerIterator<'a> {
    sequence: &'a [u8],
    m: usize,
    position: usize,
}

impl<'a> MmerIterator<'a> {
    fn new(sequence: &'a [u8], m: usize) -> Self {
        Self {
            sequence,
            m,
            position: 0,
        }
    }
}

impl<'a> Iterator for MmerIterator<'a> {
    type Item = MinimizerInfos<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position + self.m <= self.sequence.len() {
            let mmer_seq = &self.sequence[self.position..self.position + self.m];
            let mmer = OrientedSequence::new(mmer_seq, is_canonical(mmer_seq));
            let info = MinimizerInfos {
                mmer,
                position: self.position,
            };
            self.position += 1;
            Some(info)
        } else {
            None
        }
    }
}

struct MinimizerIterator<'a> {
    mmer_iter: MmerIterator<'a>,
    previous_was_none: bool,
    queue: MinimizerQueue<MinimizerInfos<'a>>,
}

impl<'a> MinimizerIterator<'a> {
    fn new(mut mmer_iter: MmerIterator<'a>, k: usize, m: usize) -> Self {
        let w = k - m + 1;
        let previous_was_none = false;
        let mut queue = MinimizerQueue::new(w);
        for mmer in mmer_iter.by_ref().take(w) {
            queue.insert(mmer);
        }
        Self {
            mmer_iter,
            previous_was_none,
            queue,
        }
    }
}

impl<'a> Iterator for MinimizerIterator<'a> {
    type Item = MinimizerInfos<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(mmer_info) = self.mmer_iter.next() {
            let minimizer = self.queue.get_min();
            self.queue.insert(mmer_info);
            Some(minimizer)
        } else if !self.previous_was_none {
            self.previous_was_none = true;
            Some(self.queue.get_min())
        } else {
            None
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

pub struct SuperkmerIterator<'a> {
    sequence: &'a [u8],
    minimizer_iter: std::iter::Peekable<MinimizerIterator<'a>>,
    k: usize,
    m: usize,
    current_minimizer: Option<MinimizerInfos<'a>>,
    start_of_superkmer_as_read: usize,
    end_of_superkmer_as_read: usize,
    i: usize,
}

impl<'a> SuperkmerIterator<'a> {
    fn new(sequence: &'a [u8], minimizer_iter: MinimizerIterator<'a>, k: usize, m: usize) -> Self {
        let mut minimizer_iter = minimizer_iter.peekable();
        let current_minimizer = minimizer_iter.next();

        Self {
            sequence,
            minimizer_iter,
            k,
            m,
            current_minimizer,
            start_of_superkmer_as_read: 0,
            end_of_superkmer_as_read: k,
            i: 1,
        }
    }
}

impl<'a> Iterator for SuperkmerIterator<'a> {
    type Item = Superkmer<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let current_minimizer = match self.current_minimizer {
            Some(ref minimizer) => minimizer,
            None => return None,
        };
        // TODO discuss style
        // let Some(current_minimizer) = self.current_minimizer else { return None };

        for candidate_minimizer in self.minimizer_iter.by_ref() {
            self.i += 1;

            if current_minimizer.position == candidate_minimizer.position {
                self.end_of_superkmer_as_read += 1;
            } else {
                let superkmer = Superkmer::new(
                    self.sequence,
                    current_minimizer.position,
                    current_minimizer.position + self.m,
                    self.start_of_superkmer_as_read,
                    self.end_of_superkmer_as_read,
                    current_minimizer.mmer.is_same_orientation(),
                );

                self.start_of_superkmer_as_read = self.i - 1;
                self.end_of_superkmer_as_read = self.i - 1 + self.k;
                self.current_minimizer = Some(candidate_minimizer);

                return Some(superkmer);
            }
        }

        if self.start_of_superkmer_as_read < self.sequence.len() {
            let superkmer = Superkmer::new(
                self.sequence,
                current_minimizer.position,
                current_minimizer.position + self.m,
                self.start_of_superkmer_as_read,
                self.end_of_superkmer_as_read,
                current_minimizer.mmer.is_same_orientation(),
            );
            self.start_of_superkmer_as_read = self.sequence.len(); // ensure no more superkmers are returned
            return Some(superkmer);
        }

        None
    }
}

pub fn compute_superkmers_linear_streaming(
    sequence: &[u8],
    k: usize,
    m: usize,
) -> Option<SuperkmerIterator> {
    if sequence.len() < k {
        None
    } else {
        let mmers = MmerIterator::new(sequence, m);
        let minimizer_iter = MinimizerIterator::new(mmers, k, m);
        // let minimizers: Vec<MinimizerInfos> = minimizer_iter.collect();
        let superkmer_iter = SuperkmerIterator::new(sequence, minimizer_iter, k, m);
        Some(superkmer_iter)
    }
}

// // Function to compute superkmers
// pub fn compute_superkmers_linear(sequence: &str, k: usize, m: usize) -> Vec<Superkmer> {
//     let superkmers = Vec::with_capacity(sequence.len());
//     if sequence.len() < k {
//         return superkmers;
//     }

//     let mmers = MmerIterator::new(sequence, m);
//     let minimizer_iter = MinimizerIterator::new(mmers, k, m);
//     // let minimizers: Vec<MinimizerInfos> = minimizer_iter.collect();
//     let superkmer_iter = SuperkmerIterator::new(sequence, minimizer_iter, k, m);
//     let superkmers = superkmer_iter.collect();

//     superkmers
// }

#[cfg(test)]
mod tests {
    // use crate::superkmer;

    use super::*;

    #[test]
    fn test_reverse_complement() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTGCA".as_bytes()).collect();
        assert_eq!(revcomp, b"TGCACTGCACAGT");
    }

    #[test]
    fn test_reverse_complement_n() {
        let revcomp: Vec<u8> = reverse_complement("ACTGTGCAGTNNGNCA".as_bytes()).collect();
        assert_eq!(revcomp, b"TGNCNNACTGCACAGT");
    }

    #[test]
    fn test_get_canonical_kmer() {
        let kmer = "ACTGCGATGACGCAGATAGCAGATAGC".as_bytes();
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
        let kmer = "GCTATCTGCTATCTGCGTCATCGCAGT".as_bytes();
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
        let seq = "AGCAGCTAGCATTTTTGCAGT".as_bytes();
        let superkmers: Vec<Superkmer> = compute_superkmers_linear_streaming(seq, 16, 5)
            .unwrap()
            .collect();
        assert_eq!(
            // ACTGCAAAAATGCTAGCTGCT
            superkmers,
            vec![Superkmer::new(seq, 11, 11 + 5, 0, seq.len(), false)]
        );
    }
    // encoding and decoding
    #[test]
    fn test_compute_superkmers2() {
        // same minimizer further away (AAAAA)
        let seq = "AGCAGCTAGCATTTTTGCAGAAAAACC".as_bytes();
        let superkmers: Vec<Superkmer> = compute_superkmers_linear_streaming(seq, 16, 5)
            .unwrap()
            .collect();
        assert_eq!(
            superkmers,
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
}
