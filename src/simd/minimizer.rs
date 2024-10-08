use crate::ringbuf::RingBuf;
use crate::simd::linearize_backup;
use crate::simd::nthash::{nthash32_par_it, nthash32_scalar_it};
use crate::simd::packed::IntoBpIterator;
use core::array::from_fn;
use core::cmp::Ordering;
use itertools::Itertools;
use std::collections::VecDeque;
use wide::u32x8 as S;

const POS_BITS: u32 = 13;
const MAX_POS: u32 = (1 << POS_BITS) - 1;
const POS_MASK: u32 = (1 << POS_BITS) - 1;
const VAL_MASK: u32 = u32::MAX ^ POS_MASK;

#[inline(always)]
fn sliding_min_tie_scalar_it(
    it: impl ExactSizeIterator<Item = u32>,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    assert!(w > 0);
    assert!(w < (1 << 15), "This method is not tested for large w.");
    assert!(it.size_hint().0 < (1 << 32));
    let max_pos = MAX_POS;
    let pos_mask = POS_MASK;
    let val_mask = VAL_MASK;
    let mut prefix_min = u32::MAX;
    let mut prefix_min_masked = val_mask;
    let mut ring_buf = RingBuf::new(w, prefix_min);
    let mut tie_ring_buf = RingBuf::new(w, 0);
    let mut pos = 0;
    let mut pos_offset = 0;
    let bool_u32 = {
        #[inline(always)]
        |x: bool| {
            if x {
                u32::MAX
            } else {
                0
            }
        }
    };

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                let delta = ((1 << POS_BITS) - 2 - w) as u32;
                pos -= delta;
                prefix_min -= delta;
                pos_offset += delta;
                for x in &mut *ring_buf {
                    *x -= delta;
                }
            }
            let elem_masked = val << POS_BITS;
            let elem = elem_masked | pos;
            ring_buf.push(elem);
            let keep_min = prefix_min <= elem;
            let prefix_tie = (bool_u32(keep_min) & tie_ring_buf.peek_prev())
                | bool_u32(elem_masked == prefix_min_masked);
            tie_ring_buf.push(prefix_tie);
            prefix_min = prefix_min.min(elem);
            prefix_min_masked = prefix_min & val_mask;
            pos += 1;
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                let mut suffix_min = ring_buf[w - 1];
                tie_ring_buf[w - 1] = 0;
                for i in (0..w - 1).rev() {
                    let keep_min = suffix_min <= ring_buf[i];
                    tie_ring_buf[i] = (bool_u32(keep_min) & tie_ring_buf[i + 1])
                        | bool_u32((ring_buf[i] & val_mask) == (suffix_min & val_mask));
                    suffix_min = suffix_min.min(ring_buf[i]);
                    ring_buf[i] = suffix_min;
                }
                prefix_min = u32::MAX;
                prefix_min_masked = val_mask;
            }
            let suffix_min = *ring_buf.peek();
            let suffix_min_masked = suffix_min & val_mask;
            let suffix_tie = *tie_ring_buf.peek();
            let tie = bool_u32(prefix_min_masked == suffix_min_masked)
                | bool_u32(prefix_min_masked < suffix_min_masked) & prefix_tie
                | bool_u32(prefix_min_masked > suffix_min_masked) & suffix_tie;
            ((prefix_min.min(suffix_min) & pos_mask) + pos_offset) | tie
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

#[inline(always)]
fn sliding_min_tie_par_it(
    it: impl ExactSizeIterator<Item = S>,
    w: usize,
) -> impl ExactSizeIterator<Item = S> {
    assert!(w > 0);
    assert!(w < (1 << 15), "This method is not tested for large w.");
    assert!(it.size_hint().0 * 8 < (1 << 32));
    let max_pos = S::splat(MAX_POS);
    let pos_mask = S::splat(POS_MASK);
    let val_mask = S::splat(VAL_MASK);
    let mut prefix_min = S::MAX;
    let mut prefix_min_masked = val_mask;
    let mut ring_buf = RingBuf::new(w, prefix_min);
    let mut tie_ring_buf = RingBuf::new(w, S::ZERO);
    let mut pos = S::ZERO;
    let mut pos_offset: S =
        from_fn(|l| (l * (it.size_hint().0.saturating_sub(w - 1))) as u32).into();

    let mut it = it.map(
        #[inline(always)]
        move |val| {
            // Make sure the position does not interfere with the hash value.
            if pos == max_pos {
                let delta = S::splat((1 << POS_BITS) - 2 - w as u32);
                pos -= delta;
                prefix_min -= delta;
                pos_offset += delta;
                for x in &mut *ring_buf {
                    *x -= delta;
                }
            }
            let elem_masked = val << POS_BITS;
            let elem = elem_masked | pos;
            ring_buf.push(elem);
            let keep_min = !prefix_min.cmp_gt(elem);
            let prefix_tie =
                (keep_min & tie_ring_buf.peek_prev()) | elem_masked.cmp_eq(prefix_min_masked);
            tie_ring_buf.push(prefix_tie);
            prefix_min = prefix_min.min(elem);
            prefix_min_masked = prefix_min & val_mask;
            pos += S::ONE;
            // After a chunk has been filled, compute suffix minima.
            if ring_buf.idx() == 0 {
                let mut suffix_min = ring_buf[w - 1];
                tie_ring_buf[w - 1] = S::ZERO;
                for i in (0..w - 1).rev() {
                    let keep_min = !suffix_min.cmp_gt(ring_buf[i]);
                    tie_ring_buf[i] = (keep_min & tie_ring_buf[i + 1])
                        | (ring_buf[i] & val_mask).cmp_eq(suffix_min & val_mask);
                    suffix_min = suffix_min.min(ring_buf[i]);
                    ring_buf[i] = suffix_min;
                }
                prefix_min = S::MAX;
                prefix_min_masked = val_mask;
            }
            let suffix_min = *ring_buf.peek();
            let suffix_min_masked = suffix_min & val_mask;
            let suffix_tie = *tie_ring_buf.peek();
            let tie = prefix_min_masked.cmp_eq(suffix_min_masked)
                | prefix_min_masked.cmp_lt(suffix_min_masked) & prefix_tie
                | prefix_min_masked.cmp_gt(suffix_min_masked) & suffix_tie;
            ((prefix_min.min(suffix_min) & pos_mask) + pos_offset) | tie
        },
    );
    // This optimizes better than it.skip(w-1).
    it.by_ref().take(w - 1).for_each(drop);
    it
}

/// Returns the minimizer of a window using a naive linear scan.
/// Uses NT hash with canonical hashes when `RC` is true.
pub fn minimizer_window_naive<const RC: bool>(seq: impl IntoBpIterator, k: usize) -> usize {
    nthash32_scalar_it::<RC>(seq, k)
        .map(|x| x & VAL_MASK)
        .position_min()
        .unwrap()
}

/// Returns an iterator over the absolute positions of the minimizers of a sequence.
/// Returns one value for each window of size `w+k-1` in the input. Use
/// `Itertools::dedup()` to obtain the distinct positions of the minimizers.
///
/// Prefer `minimizer_simd_it` that internally used SIMD, or `minimizer_par_it` if it works for you.
pub fn minimizer_scalar_it<const RC: bool>(
    seq: impl IntoBpIterator,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    let it = nthash32_scalar_it::<RC>(seq, k);
    sliding_min_tie_scalar_it(it, w)
}

fn backup_minimizer_scalar_it<const RC: bool>(
    seq: impl IntoBpIterator,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    let mut it = nthash32_scalar_it::<RC>(seq, k)
        .map(|val| val << POS_BITS)
        .enumerate();
    let mut min_val = it.next().unwrap().1;
    let mut min_pos = VecDeque::with_capacity(w);
    min_pos.push_back(0);
    let w = w as u32;
    for i in 1..(w - 1) {
        let val = it.next().unwrap().1;
        match val.cmp(&min_val) {
            Ordering::Less => {
                min_val = val;
                min_pos.clear();
                min_pos.push_back(i);
            }
            Ordering::Equal => min_pos.push_back(i),
            Ordering::Greater => (),
        }
    }
    let mut min_center_idx = 0;
    let mut min_center_pos = min_pos[0];
    for (j, &pos) in min_pos.iter().enumerate() {
        let left_dist = min_center_pos;
        let right_dist = w - 1 - pos;
        match left_dist.cmp(&right_dist) {
            Ordering::Less => {
                min_center_idx = j;
                min_center_pos = pos;
            }
            Ordering::Equal => {
                let canonical = seq.get_bp((w - 1 - w / 2) as usize) & 0b10 == 0;
                if canonical {
                    min_center_idx = j;
                    min_center_pos = pos;
                }
                break;
            }
            Ordering::Greater => break,
        }
    }
    it.map(move |(i, val)| {
        if i as u32 - w == min_pos[0] {
            min_pos.pop_front();
            min_center_idx = min_center_idx.saturating_sub(1);
        }
        if val == min_val {
            min_pos.push_back(i as u32);
        }
        if min_center_idx + 1 < min_pos.len() {
            let pos = min_pos[min_center_idx + 1];
            let left_dist = min_center_pos + w - 1 - i as u32;
            let right_dist = i as u32 - pos;
            match left_dist.cmp(&right_dist) {
                Ordering::Less => {
                    min_center_idx += 1;
                    min_center_pos = pos;
                }
                Ordering::Equal => {
                    let canonical = seq.get_bp((i as u32 - w / 2) as usize) & 0b10 == 0;
                    if canonical {
                        min_center_idx += 1;
                        min_center_pos = pos;
                    }
                }
                Ordering::Greater => (),
            }
        }
        min_center_pos
    })
}

/// Returns an iterator over the absolute positions of the minimizers of a sequence.
/// Returns one value for each window of size `w+k-1` in the input. Use
/// `Itertools::dedup()` to obtain the distinct positions of the minimizers.
///
/// This splits the windows of the sequence into chunks of 2^13.
/// Minimizers for each chunk are eagerly computed using 8 parallel streams using SIMD using `minimizers_par_it`.
/// Then returns a linear iterator over the buffer.
/// Once the buffer runs out, the next chunk is computed.
///
/// NOTE: This method is ~4x slower than the minimizer computation itself, and
///       only ~2x faster than the scalar version. Mostly because shuffling memory is slow.
/// TODO: Fix this.
pub fn minimizer_simd_it<const RC: bool>(
    seq: impl IntoBpIterator,
    k: usize,
    w: usize,
) -> impl ExactSizeIterator<Item = u32> {
    linearize_backup::linearize_with_offset(
        seq,
        k + w - 1,
        move |seq| minimizer_par_it::<RC>(seq, k, w),
        move |seq| backup_minimizer_scalar_it::<RC>(seq, k, w),
    )
}

/// Split the windows of the sequence into 8 chunks of equal length ~len/8.
/// Then return the positions of the minimizers of each of them in parallel using SIMD,
/// and return the remaining few using the second iterator.
// TODO: Take a hash function as argument.
pub fn minimizer_par_it<const RC: bool>(
    seq: impl IntoBpIterator,
    k: usize,
    w: usize,
) -> (
    impl ExactSizeIterator<Item = S>,
    impl ExactSizeIterator<Item = u32>,
) {
    let (par_head, tail) = nthash32_par_it::<RC>(seq, k, w);
    let par_head = sliding_min_tie_par_it(par_head, w);
    let offset = 8 * par_head.size_hint().0 as u32;
    let tail = sliding_min_tie_scalar_it(tail, w).map(move |pos| pos.saturating_add(offset));
    (par_head, tail)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::simd::packed::Packed;
    use rand::random;
    use std::{cell::LazyCell, iter::once};

    const BYTE_SEQ: LazyCell<Vec<u8>> =
        LazyCell::new(|| (0..1024 * 1024).map(|_| random::<u8>() & 0b110).collect());
    const PACKED_SEQ: LazyCell<Vec<u8>> =
        LazyCell::new(|| (0..1024 * 1024 / 4).map(|_| random::<u8>()).collect());

    #[test]
    fn test_sliding_min_tie_scalar_1() {
        let w = 2;
        let val = [2, 1, 2, 2, 1, 1, 3, 1, 3];
        let scalar = sliding_min_tie_scalar_it(val.iter().copied(), w).collect::<Vec<_>>();
        assert_eq!(scalar, [1, 1, u32::MAX, 4, u32::MAX, 5, 7, 7]);
    }

    #[test]
    fn test_sliding_min_tie_scalar_2() {
        let w = 2;
        let val = [2, 2, 2, 0, 1, 0, 3, 0, 2];
        let scalar = sliding_min_tie_scalar_it(val.iter().copied(), w).collect::<Vec<_>>();
        assert_eq!(scalar, [u32::MAX, u32::MAX, 3, 3, 5, 5, 7, 7]);
    }

    #[test]
    fn test_sliding_min_tie_scalar_3() {
        let w = 3;
        let val = [3, 0, 0, 3, 3, 3, 0, 3, 0, 2];
        let scalar = sliding_min_tie_scalar_it(val.iter().copied(), w).collect::<Vec<_>>();
        assert_eq!(scalar, [u32::MAX, u32::MAX, 2, u32::MAX, 6, 6, u32::MAX, 8]);
    }

    #[test]
    fn test_sliding_min_tie_scalar_4() {
        let w = 3;
        let val = [0, 1, 1, 2, 0, 2, 3, 3, 0, 2];
        let scalar = sliding_min_tie_scalar_it(val.iter().copied(), w).collect::<Vec<_>>();
        assert_eq!(scalar, [0, u32::MAX, 4, 4, 4, 5, 8, 8]);
    }

    #[test]
    fn test_sliding_min_tie_par_1() {
        let k = 1;
        let w = 3;
        let val = [0, 1, 1, 2, 0, 2, 3, 3, 0, 2]
            .map(|x| ((x & 1) << 1) | ((x ^ 2) >> 1)) // 1 3 0 2
            .map(|x| x << 1);
        let (par_head, tail) = minimizer_par_it::<false>(val.as_slice(), k, w);
        let par_head = par_head.collect::<Vec<_>>();
        let parallel_iter = (0..8)
            .flat_map(|l| par_head.iter().map(move |x| x.as_array_ref()[l]))
            .chain(tail)
            .collect::<Vec<_>>();
        assert_eq!(parallel_iter, [0, u32::MAX, 4, 4, 4, 5, 8, 8]);
    }

    #[test]
    fn parallel_iter_byte() {
        let seq = &**BYTE_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128)) {
                    let seq = seq.sub_slice(0, len);
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    let (par_head, tail) = minimizer_par_it::<false>(seq, k, w);
                    let par_head = par_head.collect::<Vec<_>>();
                    let parallel_iter = (0..8)
                        .flat_map(|l| par_head.iter().map(move |x| x.as_array_ref()[l]))
                        .chain(tail)
                        .collect::<Vec<_>>();
                    assert_eq!(scalar, parallel_iter, "k={k}, w={w}, len={len}");
                }
            }
        }
    }

    #[test]
    fn parallel_iter_packed() {
        let seq = Packed {
            seq: &*PACKED_SEQ,
            offset: 0,
            len: 1024 * 1024,
        };
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128)) {
                    let seq = seq.sub_slice(0, len);
                    let scalar = minimizer_scalar_it::<false>(seq, k, w).collect::<Vec<_>>();
                    let (par_head, tail) = minimizer_par_it::<false>(seq, k, w);
                    let par_head = par_head.collect::<Vec<_>>();
                    let parallel_iter = (0..8)
                        .flat_map(|l| par_head.iter().map(move |x| x.as_array_ref()[l]))
                        .chain(tail)
                        .collect::<Vec<_>>();
                    assert_eq!(scalar, parallel_iter, "k={k}, w={w}, len={len}");
                }
            }
        }
    }

    #[test]
    fn linearized_increasing() {
        let seq = &**BYTE_SEQ;
        for k in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
            for w in [1, 2, 3, 4, 5, 31, 32, 33, 63, 64, 65] {
                for len in (0..100).chain(once(1024 * 128 + 765)) {
                    let seq = seq.sub_slice(0, len);
                    let mut prev = 0;
                    let simd = minimizer_simd_it::<false>(seq, k, w).collect::<Vec<_>>();
                    for &pos in simd.iter() {
                        assert!(prev <= pos, "k={k}, w={w}, len={len}");
                        prev = pos;
                    }
                }
            }
        }
    }
}
