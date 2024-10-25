use std::alloc::{alloc_zeroed, dealloc, Layout};
use std::rc::Rc;
use std::sync::RwLock;

use crate::subsequence::{NoBitPacked, Subsequence};

/// Wrapper around a pointer.
/// YOU HAVE TO DEALLOCATE IT BEFORE DROPPING IT
pub struct ArrayOfHyperkmer {
    ptr: Rc<RwLock<*mut u64>>,

    #[cfg(any(debug_assertions, test))]
    is_dealloc: bool,
}

// TODO is it safe ?
unsafe impl Send for ArrayOfHyperkmer {}
unsafe impl Sync for ArrayOfHyperkmer {}

impl ArrayOfHyperkmer {
    /// Creates a new `Self` of a given `size`.
    /// `size` should be the number of `u64` fitting into `Self`.
    pub fn new(size: usize) -> Self {
        // Create a memory layout for the array
        let layout = Layout::array::<u64>(size).expect("Failed to create layout");

        // Allocate the memory
        let ptr = unsafe { alloc_zeroed(layout) };

        // Check for null pointer
        if ptr.is_null() {
            panic!("Memory allocation failed");
        }

        let ptr = ptr as *mut u64;
        let ptr = Rc::new(RwLock::new(ptr));

        ArrayOfHyperkmer {
            ptr,
            #[cfg(any(debug_assertions, test))]
            is_dealloc: false,
        }
    }

    /// Constructs a `Self` from a u64 iterator.
    pub fn from_u64_iter<I>(iter: I, size: usize) -> Self
    where
        I: IntoIterator<Item = u64>,
    {
        let simple_vec = Self::new(size);

        // Copy elements from iterator into the allocated memory
        for (i, value) in iter.into_iter().enumerate() {
            if i >= size {
                panic!("Iterator provided more elements than expected size");
            }
            // DEBUG: no need to use `write` as this is a new vector ?
            let ptr = simple_vec.ptr.read().expect("could not acquire read lock");
            unsafe {
                *ptr.add(i) = value;
            }
        }

        simple_vec
    }

    /// Computes a slice.
    /// `len` is the number of `u64` in the slice.
    pub fn as_slice(&self, len: usize) -> &[u64] {
        let ptr = self.ptr.read().expect("could not acquire read lock");
        unsafe { std::slice::from_raw_parts(*ptr, len * 64) }
    }

    pub fn dump(&mut self, len: usize, start: usize, end: usize, subseq: Subsequence<NoBitPacked>) {
        let ptr = self.ptr.read().expect("could not acquire read lock");
        let slice: &mut [u64] = unsafe { std::slice::from_raw_parts_mut(*ptr, len * 64) };
        subseq.dump_as_2bits(&mut slice[start..end]);
    }

    // DEBUG this lets the memory be accessed out of the lock (but not modified)
    // How to prevent memory from being modified ?
    pub fn as_u64_slice(&self, len: usize) -> &[u64] {
        let ptr = self.ptr.read().expect("could not acquire read lock");
        unsafe { std::slice::from_raw_parts(*ptr, len) }
    }

    pub fn dealloc(&mut self, size: usize) {
        // Create a memory layout for the array
        let layout = Layout::array::<u8>(size).expect("Failed to create layout");

        // Deallocate the memory
        let ptr = self.ptr.read().expect("could not acquire write lock");
        unsafe {
            dealloc(*ptr as *mut u8, layout);
        }

        #[cfg(any(debug_assertions, test))]
        {
            self.is_dealloc = true;
        }
    }
}

impl Drop for ArrayOfHyperkmer {
    fn drop(&mut self) {
        #[cfg(any(debug_assertions, test))]
        {
            if !self.is_dealloc {
                panic!()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn forgot_alloc() {
        ArrayOfHyperkmer::new(100);
    }
}
