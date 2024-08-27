use std::alloc::{alloc_zeroed, dealloc, Layout};
use std::rc::Rc;
use std::sync::RwLock;

use crate::superkmer::{NoBitPacked, SubsequenceMetadata};

/// Wrapper around a pointer.
/// YOU HAVE TO DEALLOCATE IT BEFORE DROPPING IT
pub struct SimpleVec {
    ptr: Rc<RwLock<*mut u64>>,

    #[cfg(debug_assertions)]
    is_dealloc: bool,
}

unsafe impl Send for SimpleVec {}
unsafe impl Sync for SimpleVec {}

impl SimpleVec {
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

        SimpleVec {
            ptr,
            #[cfg(debug_assertions)]
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

    // TODO this seems to work, but is this UB ?
    pub fn as_slice(&self, len: usize) -> &[u8] {
        let ptr = self.ptr.read().expect("could not acquire read lock");
        unsafe { std::slice::from_raw_parts(*ptr as *mut u8, len * 8) }
    }

    pub fn dump(
        &mut self,
        len: usize,
        start: usize,
        end: usize,
        subseq: SubsequenceMetadata<NoBitPacked>,
    ) {
        let ptr = self.ptr.write().expect("could not acquire write lock");
        let slice: &mut [u8] = unsafe { std::slice::from_raw_parts_mut(*ptr as *mut u8, len * 8) };
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
        let ptr = self.ptr.write().expect("could not acquire write lock");
        unsafe {
            dealloc(*ptr as *mut u8, layout);
        }

        #[cfg(debug_assertions)]
        {
            self.is_dealloc = true;
        }
    }
}

impl Drop for SimpleVec {
    fn drop(&mut self) {
        #[cfg(debug_assertions)]
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
        SimpleVec::new(100);
    }
}
