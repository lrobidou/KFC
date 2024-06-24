use std::alloc::{alloc, dealloc, Layout};
use std::ptr;

/// Wrapper around a pointer.
/// YOU HAVE TO DEALLOCATE IT BEFORE DROPPING IT
pub struct SimpleVec {
    ptr: *mut u8,

    #[cfg(debug_assertions)]
    is_dealloc: bool,
}

impl SimpleVec {
    pub fn new(size: usize) -> Self {
        // Create a memory layout for the array
        let layout = Layout::array::<u8>(size).expect("Failed to create layout");

        // Allocate the memory
        let ptr = unsafe { alloc(layout) };

        // Check for null pointer
        if ptr.is_null() {
            panic!("Memory allocation failed");
        }

        // Initialize the memory with zeros
        unsafe {
            ptr::write_bytes(ptr, 0, size);
        }

        SimpleVec {
            ptr,
            #[cfg(debug_assertions)]
            is_dealloc: false,
        }
    }

    #[cfg(test)]
    pub fn from_iter<I>(iter: I, size: usize) -> Self
    where
        I: IntoIterator<Item = u8>,
    {
        let simple_vec = Self::new(size);

        // Copy elements from iterator into the allocated memory
        for (i, byte) in iter.into_iter().enumerate() {
            if i >= size {
                panic!("Iterator provided more elements than expected size");
            }
            unsafe {
                *simple_vec.ptr.add(i) = byte;
            }
        }

        simple_vec
    }

    pub fn as_slice(&self, len: usize) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr, len) }
    }

    pub fn as_mut_slice(&mut self, len: usize) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr, len) }
    }

    pub fn dealloc(&mut self, size: usize) {
        // Create a memory layout for the array
        let layout = Layout::array::<u8>(size).expect("Failed to create layout");

        // Deallocate the memory
        unsafe {
            dealloc(self.ptr, layout);
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

    #[test]
    fn test_from_iter() {
        let kmer = "adfygkhbalzejchv";
        let k = kmer.len();
        let mut sv = SimpleVec::from_iter(kmer.as_bytes().iter().copied(), k);
        assert_eq!(kmer, &String::from_utf8(sv.as_slice(k).into()).unwrap());
        sv.dealloc(k);
    }
}
