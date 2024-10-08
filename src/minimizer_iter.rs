use crate::simd::minimizer::minimizer_simd_it;
use core::cmp::min;
use core::hash::Hash;
use itertools::Itertools;
use minimizer_queue::MinimizerQueue;
use num_traits::{AsPrimitive, PrimInt};
use std::collections::VecDeque;

pub struct MinimizerItem {
    pub position: usize,
    pub window_start: usize,
}

pub fn simd_minimizer_iter(
    seq: &[u8],
    minimizer_size: usize,
    width: usize,
) -> impl Iterator<Item = MinimizerItem> + '_ {
    assert_eq!(
        width % 2,
        1,
        "width must be odd to break ties between multiple minimizers"
    );
    let it = minimizer_simd_it::<true>(seq, minimizer_size, width);
    it.enumerate()
        .dedup_by(|(_, p), (_, q)| p == q)
        .map(|(window_start, pos)| MinimizerItem {
            position: pos as usize,
            window_start,
        })
}

/// An iterator over the canonical minimizers of a sequence and their positions
/// with the starting positions of the associated windows and a boolean indicating a reverse complement.
/// It requires an odd width to break ties between multiple minimizers.
pub struct CanonicalMinimizerIterator<'a, T: PrimInt + Hash = u64> {
    pub(crate) seq: &'a [u8],
    pub(crate) queue: MinimizerQueue<T>,
    pub(crate) width: usize,
    pub(crate) mmer: T,
    pub(crate) rc_mmer: T,
    pub(crate) mmer_mask: T,
    pub(crate) rc_mmer_shift: usize,
    pub(crate) is_rc: VecDeque<bool>,
    pub(crate) encoding: [u8; 256],
    pub(crate) rc_encoding: [u8; 256],
    pub(crate) base_width: usize,
    pub(crate) min_pos: (T, usize, usize, bool),
    pub(crate) end: usize,
}

impl<'a, T: PrimInt + Hash> CanonicalMinimizerIterator<'a, T> {
    #[deprecated(note = "This is a deprecated implementation, the SIMD version is recommended.")]
    pub fn new(seq: &'a [u8], minimizer_size: usize, width: u16) -> Self {
        let queue = MinimizerQueue::new(width);
        let width = width as usize;
        assert_eq!(
            width % 2,
            1,
            "width must be odd to break ties between multiple minimizers"
        );
        let mut encoding = [0u8; 256];
        encoding[b'A' as usize] = 0b00;
        encoding[b'a' as usize] = 0b00;
        encoding[b'C' as usize] = 0b01;
        encoding[b'c' as usize] = 0b01;
        encoding[b'G' as usize] = 0b10;
        encoding[b'g' as usize] = 0b10;
        encoding[b'T' as usize] = 0b11;
        encoding[b't' as usize] = 0b11;
        let mut rc_encoding = encoding;
        rc_encoding.swap(b'A' as usize, b'T' as usize);
        rc_encoding.swap(b'a' as usize, b't' as usize);
        rc_encoding.swap(b'C' as usize, b'G' as usize);
        rc_encoding.swap(b'c' as usize, b'g' as usize);

        encoding[b'N' as usize] = 0b00;
        encoding[b'n' as usize] = 0b00;
        rc_encoding[b'N' as usize] = 0b11;
        rc_encoding[b'n' as usize] = 0b11;
        Self {
            seq,
            queue,
            width,
            mmer: T::zero(),
            rc_mmer: T::zero(),
            mmer_mask: (T::one() << (2 * minimizer_size)) - T::one(),
            rc_mmer_shift: 2 * (minimizer_size - 1),
            is_rc: VecDeque::with_capacity(width),
            encoding,
            rc_encoding,
            base_width: width + minimizer_size - 1,
            end: width + minimizer_size - 1,
            min_pos: (T::zero(), 0, 0, false),
        }
    }

    #[inline]
    fn window_not_canonical(&self) -> bool {
        self.is_rc[self.width / 2]
    }
}

impl<'a, T: PrimInt + Hash + 'static> Iterator for CanonicalMinimizerIterator<'a, T>
where
    u8: AsPrimitive<T>,
{
    type Item = (T, usize, usize, bool);

    fn next(&mut self) -> Option<Self::Item> {
        if self.queue.is_empty() {
            if self.base_width > self.seq.len() {
                return None;
            }
            for i in 0..(self.base_width - self.width) {
                self.mmer = (self.mmer << 2)
                    | (unsafe { self.encoding.get_unchecked(self.seq[i] as usize) }.as_());
                self.rc_mmer = (self.rc_mmer >> 2)
                    | (unsafe { self.rc_encoding.get_unchecked(self.seq[i] as usize) }.as_()
                        << self.rc_mmer_shift);
            }
            for i in (self.base_width - self.width)..self.base_width {
                self.mmer = ((self.mmer << 2) & self.mmer_mask)
                    | (unsafe { self.encoding.get_unchecked(self.seq[i] as usize) }.as_());
                self.rc_mmer = (self.rc_mmer >> 2)
                    | (unsafe { self.rc_encoding.get_unchecked(self.seq[i] as usize) }.as_()
                        << self.rc_mmer_shift);
                let canonical_mmer = min(self.mmer, self.rc_mmer);
                self.queue.insert(canonical_mmer);
                self.is_rc.push_back(canonical_mmer == self.rc_mmer);
            }
            let _min_pos = if self.queue.multiple_mins() {
                let (x, pos, tie) = self.queue.get_inner_min_pos();
                tie.map_or((x, pos), |alt| {
                    if self.window_not_canonical() {
                        alt
                    } else {
                        (x, pos)
                    }
                })
            } else {
                self.queue.get_min_pos()
            };
            self.min_pos = (_min_pos.0, _min_pos.1, 0, self.is_rc[_min_pos.1]);
        } else {
            let mut min_pos = self.min_pos;
            while self.end < self.seq.len() && min_pos.1 == self.min_pos.1 {
                self.mmer = ((self.mmer << 2) & self.mmer_mask)
                    | (unsafe { self.encoding.get_unchecked(self.seq[self.end] as usize) }.as_());
                self.rc_mmer = (self.rc_mmer >> 2)
                    | (unsafe { self.rc_encoding.get_unchecked(self.seq[self.end] as usize) }
                        .as_()
                        << self.rc_mmer_shift);
                let canonical_mmer = min(self.mmer, self.rc_mmer);
                self.queue.insert(canonical_mmer);
                self.is_rc.pop_front();
                self.is_rc.push_back(canonical_mmer == self.rc_mmer);
                self.end += 1;
                let _min_pos = if self.queue.multiple_mins() {
                    let (x, pos, tie) = self.queue.get_inner_min_pos();
                    tie.map_or((x, pos), |alt| {
                        if self.window_not_canonical() {
                            alt
                        } else {
                            (x, pos)
                        }
                    })
                } else {
                    self.queue.get_min_pos()
                };
                min_pos = (
                    _min_pos.0,
                    self.end - self.base_width + _min_pos.1,
                    self.end - self.base_width,
                    self.is_rc[_min_pos.1],
                );
            }
            if min_pos.1 == self.min_pos.1 {
                return None;
            }
            self.min_pos = min_pos;
        }
        Some(self.min_pos)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_iter_single_kmer() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACG".as_bytes();
        let k = kmer.len();
        let m = 10;
        let minimizer_iter: CanonicalMinimizerIterator<u64> =
            CanonicalMinimizerIterator::new(kmer, m, (k - m) as u16);
        let minimizer_iter = minimizer_iter.collect_vec();
        assert_eq!(minimizer_iter.len(), 1);
    }
}
