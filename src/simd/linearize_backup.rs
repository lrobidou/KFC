// TODO: The linearized functions are quite slow.

use super::packed::IntoBpIterator;
use super::packed::{L, S};
use core::array::from_fn;
use core::cmp::Ordering;
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;
use core::mem::MaybeUninit;

const CHUNK_SIZE: usize = 1 << 12;

pub fn linearize_with_offset<
    BI: IntoBpIterator,
    I1: ExactSizeIterator<Item = S>,
    I2: ExactSizeIterator<Item = u32>,
    I3: ExactSizeIterator<Item = u32>,
>(
    seq: BI,
    context: usize,
    par_it: impl Fn(BI) -> (I1, I2),
    backup_it: impl Fn(BI) -> I3,
) -> impl ExactSizeIterator<Item = u32> {
    Linear::<true, _, _, _, _, _, _>::new(seq, context, par_it, backup_it)
}

struct Linear<
    const OFFSET: bool,
    BI: IntoBpIterator,
    I1: ExactSizeIterator<Item = S>,
    I2: ExactSizeIterator<Item = u32>,
    I3: ExactSizeIterator<Item = u32>,
    ParIt: Fn(BI) -> (I1, I2),
    BackupIt: Fn(BI) -> I3,
> {
    seq: BI,
    par_it: ParIt,
    backup_it: BackupIt,
    context: usize,

    chunk_size: usize,
    par_size: usize,
    num_chunks: usize,

    cache: Vec<[u32; L]>,
    tail_cache: Vec<u32>,

    /// current chunk
    c: usize,
    /// current lane in chunk
    l: usize,
    /// current index
    i: usize,

    /// number of lanes in the current chunk.
    max_l: usize,
    /// length of the current lane
    max_i: usize,
}

impl<
        const OFFSET: bool,
        BI: IntoBpIterator,
        I1: ExactSizeIterator<Item = S>,
        I2: ExactSizeIterator<Item = u32>,
        I3: ExactSizeIterator<Item = u32>,
        ParIt: Fn(BI) -> (I1, I2),
        BackupIt: Fn(BI) -> I3,
    > Linear<OFFSET, BI, I1, I2, I3, ParIt, BackupIt>
{
    fn new(
        seq: BI,
        context: usize,
        par_it: ParIt,
        backup_it: BackupIt,
    ) -> Linear<OFFSET, BI, I1, I2, I3, ParIt, BackupIt> {
        if OFFSET {
            assert!(
                seq.len() < u32::MAX as usize,
                "sequence of length {} is too long for u32 indices",
                seq.len()
            );
        }
        assert!(context > 0);
        let len = seq.len();
        let chunk_size = CHUNK_SIZE;
        assert!(chunk_size % L == 0);
        let par_size = chunk_size / L;

        let num_chunks = len.saturating_sub(context - 1).div_ceil(chunk_size);

        let cache = vec![[0; L]; par_size];
        let tail_cache = vec![];

        let mut this = Self {
            seq,
            par_it,
            backup_it,
            context,
            cache,
            tail_cache,

            chunk_size,
            par_size,
            num_chunks,

            c: usize::MAX,
            i: par_size - 1,
            l: L - 1,

            max_i: par_size,
            max_l: L,
        };
        if num_chunks == 1 {
            this.buffer_chunk(0, seq);
        }
        this
    }

    /// Fill a buffer with NT hash values by using `nthash32_par_it`.
    fn buffer_chunk(&mut self, mut offset: usize, seq: BI) {
        let (mut par_it, mut tail) = (self.par_it)(seq);
        self.cache.resize(par_it.size_hint().0, [0; L]);
        if !OFFSET {
            offset = 0;
        }
        let offset_simd = S::splat(offset as u32);
        let mut ties: [Vec<usize>; L] = from_fn(|_| Vec::new());
        for i in 0..self.cache.len() {
            unsafe {
                let pos_simd = par_it.next().unwrap();
                *self.cache.get_unchecked_mut(i) = (pos_simd + offset_simd).to_array();
                let tie_simd = pos_simd.cmp_eq(S::MAX);
                if tie_simd.any() {
                    for (j, &tie) in tie_simd.to_array().iter().enumerate() {
                        if tie == u32::MAX {
                            ties[j].push(i);
                        }
                    }
                }
            }
        }
        let offset = offset as u32;
        self.tail_cache.resize(tail.size_hint().0, 0);
        let mut tail_tie = Vec::new();
        for i in 0..self.tail_cache.len() {
            unsafe {
                let pos = tail.next().unwrap();
                *self.tail_cache.get_unchecked_mut(i) = pos + offset;
                if pos == u32::MAX {
                    tail_tie.push(i);
                }
            }
        }
        let mut tie_start = 0;
        let mut backup_it = MaybeUninit::uninit();
        for (j, tie) in ties.iter().enumerate() {
            if unlikely(!tie.is_empty()) {
                let mut prev = usize::MAX - 1;
                for &i in tie.iter() {
                    if i != prev + 1 {
                        // dbg!("tie!");
                        tie_start = j * self.cache.len() + i;
                        let len = self.cache.len() + self.context - 1;
                        backup_it =
                            MaybeUninit::new((self.backup_it)(seq.sub_slice(tie_start, len - i)));
                    }
                    self.cache[i][j] = unsafe {
                        backup_it.assume_init_mut().next().unwrap() + tie_start as u32 + offset
                    };
                    prev = i;
                }
            }
        }
        if unlikely(!tail_tie.is_empty()) {
            let mut prev = usize::MAX - 1;
            for &i in tail_tie.iter() {
                if i != prev + 1 {
                    // dbg!("tie!!");
                    tie_start = L * self.cache.len() + i;
                    let len = self.tail_cache.len() + self.context - 1;
                    backup_it =
                        MaybeUninit::new((self.backup_it)(seq.sub_slice(tie_start, len - i)));
                }
                self.tail_cache[i] = unsafe {
                    backup_it.assume_init_mut().next().unwrap() + tie_start as u32 + offset
                };
                prev = i;
            }
        }
    }

    /// Returns true when iterator is exhausted.
    #[inline(always)]
    fn next_line(&mut self) -> bool {
        self.i = 0;
        self.l += 1;
        if self.l == self.max_l {
            self.l = 0;
            self.c = self.c.wrapping_add(1);
            if self.c < self.num_chunks {
                let offset = self.c * self.chunk_size;
                self.buffer_chunk(
                    offset,
                    self.seq.sub_slice(
                        offset,
                        (self.chunk_size + self.context - 1).min(self.seq.len() - offset),
                    ),
                );
                if self.c + 1 < self.num_chunks {
                    assert_eq!(self.cache.len(), self.par_size);
                    assert!(self.tail_cache.is_empty());
                } else {
                    assert!(self.cache.len() <= self.par_size);
                    self.max_i = self.cache.len();
                    if self.max_i == 0 {
                        self.c += 1;
                    }
                }
            }
            if self.c == self.num_chunks {
                // Run one extra iteration on the tail.
                self.max_l = 1;
                self.max_i = self.tail_cache.len();
                assert!(self.tail_cache.len() <= self.par_size);
                self.cache.resize(self.max_i, [0; L]);
                for i in 0..self.max_i {
                    self.cache[i][0] = self.tail_cache[i];
                }
                if self.max_i == 0 {
                    self.c += 1;
                }
            }
            if self.c == self.num_chunks + 1 {
                return true;
            }
        }
        false
    }
}

impl<
        const OFFSET: bool,
        BI: IntoBpIterator,
        I1: ExactSizeIterator<Item = S>,
        I2: ExactSizeIterator<Item = u32>,
        I3: ExactSizeIterator<Item = u32>,
        ParIt: Fn(BI) -> (I1, I2),
        BackupIt: Fn(BI) -> I3,
    > Iterator for Linear<OFFSET, BI, I1, I2, I3, ParIt, BackupIt>
{
    type Item = u32;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.i += 1;
        if self.i == self.max_i && self.next_line() {
            return None;
        }

        Some(unsafe { *self.cache.get_unchecked(self.i).get_unchecked(self.l) })
    }

    #[inline(always)]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = match self.c.cmp(&self.num_chunks) {
            Ordering::Greater => 0,
            Ordering::Equal => self.max_i - self.i,
            Ordering::Less => {
                let total_kmers = self.seq.len() - self.context + 1;
                let processed_kmers = 1 + self.c * L * self.par_size + self.l * self.max_i + self.i;
                total_kmers - processed_kmers
            }
        };
        (len, Some(len))
    }
}

impl<
        const OFFSET: bool,
        BI: IntoBpIterator,
        I1: ExactSizeIterator<Item = S>,
        I2: ExactSizeIterator<Item = u32>,
        I3: ExactSizeIterator<Item = u32>,
        ParIt: Fn(BI) -> (I1, I2),
        BackupIt: Fn(BI) -> I3,
    > ExactSizeIterator for Linear<OFFSET, BI, I1, I2, I3, ParIt, BackupIt>
{
}
