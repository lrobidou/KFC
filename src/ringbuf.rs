#![allow(dead_code)]

/// A custom RingBuf implementation that has a fixed size `w` and wraps around.
pub struct RingBuf<V> {
    w: usize,
    idx: usize,
    data: Vec<V>,
}

impl<V: Clone> RingBuf<V> {
    #[inline(always)]
    pub fn new(w: usize, v: V) -> Self {
        assert!(w > 0);
        let data = vec![v; w];
        RingBuf { w, idx: 0, data }
    }

    /// Returns the next index to be written.
    #[inline(always)]
    pub fn idx(&self) -> usize {
        self.idx
    }

    #[inline(always)]
    pub fn prev_idx(&self) -> usize {
        if self.idx == 0 {
            self.w - 1
        } else {
            self.idx - 1
        }
    }

    #[inline(always)]
    pub fn nth_idx(&self, n: usize) -> usize {
        if self.idx + n < self.w {
            self.idx + n
        } else {
            self.idx + n - self.w
        }
    }

    #[inline(always)]
    pub fn peek(&self) -> &V {
        unsafe { self.data.get_unchecked(self.idx) }
    }

    #[inline(always)]
    pub fn peek_prev(&self) -> &V {
        unsafe { self.data.get_unchecked(self.prev_idx()) }
    }

    #[inline(always)]
    pub fn peek_nth(&self, n: usize) -> &V {
        unsafe { self.data.get_unchecked(self.nth_idx(n)) }
    }

    #[inline(always)]
    pub fn push(&mut self, v: V) {
        self.data[self.idx] = v;
        self.idx += 1;
        if self.idx == self.w {
            self.idx = 0;
        }
    }
}

/// A RingBuf can be used as a slice.
impl<V> std::ops::Deref for RingBuf<V> {
    type Target = [V];

    #[inline(always)]
    fn deref(&self) -> &[V] {
        &self.data
    }
}

/// A RingBuf can be used as a mutable slice.
impl<V> std::ops::DerefMut for RingBuf<V> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [V] {
        &mut self.data
    }
}
