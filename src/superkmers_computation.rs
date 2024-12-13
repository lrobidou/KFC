use super::superkmer::Superkmer;
use crate::minimizer_iter::CanonicalMinimizerIterator;
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

pub fn is_equal_to_its_revcomp(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };

        if xc != yc {
            return false;
        }
    }
    true
}

pub struct SuperkmerIterator<'a> {
    sequence: &'a [u8],
    minimizer_iter: std::iter::Peekable<CanonicalMinimizerIterator<'a>>,
    k: usize,
    m: usize,
    previous_minimizer: Option<(u64, usize, usize, bool)>,
}

impl<'a> SuperkmerIterator<'a> {
    fn new(
        sequence: &'a [u8],
        minimizer_iter: CanonicalMinimizerIterator<'a>,
        k: usize,
        m: usize,
    ) -> Self {
        let mut minimizer_iter = minimizer_iter.peekable();
        let previous_minimizer: Option<(u64, usize, usize, bool)> = minimizer_iter.next();
        Self {
            sequence,
            minimizer_iter,
            k,
            m,
            previous_minimizer,
        }
    }
}

impl<'a> Iterator for SuperkmerIterator<'a> {
    type Item = Superkmer<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.previous_minimizer {
            Some(ref previous_minimizer) => match self.minimizer_iter.next() {
                Some(current_minimizer) => {
                    let sk = Superkmer::new(
                        self.sequence,
                        previous_minimizer.1,
                        previous_minimizer.1 + self.m,
                        previous_minimizer.2,
                        current_minimizer.2 + self.k - 1,
                        !previous_minimizer.3,
                    );
                    self.previous_minimizer = Some(current_minimizer);
                    Some(sk)
                }
                None => {
                    let sk = Superkmer::new(
                        self.sequence,
                        previous_minimizer.1,
                        previous_minimizer.1 + self.m,
                        previous_minimizer.2,
                        self.sequence.len(),
                        !previous_minimizer.3,
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
) -> Option<SuperkmerIterator> {
    if sequence.len() < k {
        None
    } else {
        let minimizer_iter: CanonicalMinimizerIterator<u64> =
            CanonicalMinimizerIterator::new(sequence, m, (k - m + 1) as u16);
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

    #[test]
    fn test_compute_superkmers_sequence_too_short() {
        // there are no superkmers for sequence < k
        if let Some(_x) = compute_superkmers_linear_streaming("AGCAGCTAGCATTTT".as_bytes(), 16, 5) {
            panic!()
        }
    }
}
