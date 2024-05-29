// use ahash::RandomState;
// use core::hash::Hash;
use std::collections::VecDeque;

pub struct MinimizerQueue<T: Copy> {
    deq: VecDeque<(T, u8)>,
    // hash_builder: RandomState,
    w: usize,
    pos: u8,
}

impl<T: Copy + Ord> MinimizerQueue<T> {
    pub fn new(w: usize) -> Self {
        Self {
            deq: VecDeque::with_capacity(w),
            // hash_builder: RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
            pos: 0,
            w,
        }
    }

    pub fn get_min(&self) -> T {
        debug_assert!(!self.deq.is_empty(), "MinimizerQueue is empty");
        self.deq[0].0
    }

    // CAUTION: in BRR, we actually get the max
    pub fn insert(&mut self, element: T) {
        if !self.deq.is_empty() && self.deq[0].1 == self.pos {
            self.deq.pop_front();
        }
        let mut i = self.deq.len();
        while i > 0 && self.deq[i - 1].0 >= element {
            // TODO had to change this
            i -= 1;
        }
        self.deq.truncate(i);
        self.deq.push_back((element, self.pos));
        self.pos = (self.pos + 1) % self.w as u8;
    }
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;
    // use crate::brrr_kmer::{Kmer, RawKmer};

    // const K: usize = 7;
    const W: usize = 3;
    // const W: usize = K - M + 1;

    #[test]
    fn test_sliding() {
        let inputs = [1, 2, 3, 0, 7, 8, 9, 100, 3, 4, 7, 8];
        let mut queue = MinimizerQueue::<_>::new(W);
        let mut minimums: Vec<usize> = vec![];

        for (i, input) in inputs.iter().enumerate() {
            if i < W {
                queue.insert(input);
            } else {
                minimums.push(*queue.get_min());
                queue.insert(input);
            }
        }
        minimums.push(*queue.get_min());
        assert_eq!(minimums, vec![1, 0, 0, 0, 7, 8, 3, 3, 3, 4]);
    }
}
