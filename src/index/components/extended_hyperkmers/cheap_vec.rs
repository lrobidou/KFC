use std::alloc::{alloc_zeroed, dealloc, Layout};
use std::ptr;

/// Wrapper around a pointer.
/// YOU HAVE TO DEALLOCATE IT BEFORE DROPPING IT
pub struct SimpleVec {
    ptr: *mut u64,

    #[cfg(debug_assertions)]
    is_dealloc: bool,
}

impl SimpleVec {
    /// Creates a new `Self` of a give `size`.
    /// `size` should be the number a u64 fitting into `Self`.
    pub fn new(size: usize) -> Self {
        // Create a memory layout for the array
        let layout = Layout::array::<u64>(size).expect("Failed to create layout");

        // Allocate the memory
        let ptr = unsafe { alloc_zeroed(layout) };

        // Check for null pointer
        if ptr.is_null() {
            panic!("Memory allocation failed");
        }

        // Initialize the memory with zeros
        unsafe {
            ptr::write_bytes(ptr, 0, size);
        }

        SimpleVec {
            ptr: ptr as *mut u64,
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
            unsafe {
                *simple_vec.ptr.add(i) = value;
            }
        }

        simple_vec
    }

    // TODO this seems to work, but is this UB ?
    pub fn as_slice(&self, len: usize) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr as *mut u8, len * 8) }
    }

    pub fn as_mut_slice(&mut self, len: usize) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr as *mut u8, len * 8) }
    }

    pub fn as_u64_slice(&self, len: usize) -> &[u64] {
        unsafe { std::slice::from_raw_parts(self.ptr, len) }
    }

    pub fn dealloc(&mut self, size: usize) {
        // Create a memory layout for the array
        let layout = Layout::array::<u8>(size).expect("Failed to create layout");

        // Deallocate the memory
        unsafe {
            dealloc(self.ptr as *mut u8, layout);
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

    // #[test]
    // fn test_from_iter() {
    //     let kmer = "adfygkhbalzejchv";
    //     let k = kmer.len();
    //     let mut sv = SimpleVec::from_iter(vec![0, 4, 78987], k);
    //     assert_eq!(kmer, &String::from_utf8(sv.as_slice(k).into()).unwrap());
    //     sv.dealloc(k);
    // }
}