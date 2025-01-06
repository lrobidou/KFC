use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct CachedValue {
    id_bucket: usize,
    id_hk: usize,
    is_large: bool,
}

impl CachedValue {
    pub fn new(id_bucket: usize, id_hk: usize, is_large: bool) -> Self {
        Self {
            id_bucket,
            id_hk,
            is_large,
        }
    }

    pub fn get_id_bucket(&self) -> usize {
        self.id_bucket
    }

    pub fn get_id_hk(&self) -> usize {
        self.id_hk
    }

    pub fn get_is_large(&self) -> bool {
        self.is_large
    }
}

// #[repr(C)]  // TODO should I use it?
/// A cache storing a vector of references.
/// By storing a vector in the cache, the user can reuse the capacity of that vector e.g. in another part of their program,
/// even if the lifetime of the references would have prevented that.
pub struct CacheMisere<T> {
    data: *const libc::c_void, // the pointer to the vector
    len: usize,                // the length of the vector
    capacity: usize,           // the capacity of the vector
    _t: PhantomData<T>,        // keeps track of the type at compile time
}

impl<T> CacheMisere<T> {
    /// Creates a cache storing an empty vector of references.
    pub fn new() -> Self {
        Vec::new().into()
    }

    /// Internal function to be called by into
    ///
    /// # Safety
    ///
    /// Once this method is called on a `CacheMisere` instance, no other method of this particular instance should be called.
    #[allow(clippy::wrong_self_convention)]
    unsafe fn into_inner(&self) -> Vec<&T> {
        let v: Vec<&T> =
            unsafe { Vec::from_raw_parts(self.data as *mut &T, self.len, self.capacity) };
        v
    }
}

impl<T> Default for CacheMisere<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Drop for CacheMisere<T> {
    fn drop(&mut self) {
        // Safety: This is the destuctor of `CacheMisere`, so no method of this instance will be called after calling this function
        let _ = unsafe { self.into_inner() };
    }
}

impl<T> From<Vec<&T>> for CacheMisere<T> {
    /// Stores the vector in the cache, allowing it to be reused afterwards
    fn from(mut v: Vec<&T>) -> Self {
        // clear the vector from its element (does not change the capacity)
        v.clear();

        // "cast" the vector to a pointer
        let incativated = CacheMisere::<T> {
            data: v.as_ptr() as *const libc::c_void,
            len: v.len(),
            capacity: v.capacity(),
            _t: PhantomData,
        };

        // tell the compiler to not execute the destructor of the vector
        std::mem::forget(v);

        incativated
    }
}

impl<'a, T> From<CacheMisere<T>> for Vec<&'a T> {
    /// Returns the vector stored in the cache
    fn from(cache: CacheMisere<T>) -> Self {
        // Safety: cache.data was made from a Vec<&T>
        let v: Vec<&T> =
            unsafe { Vec::from_raw_parts(cache.data as *mut &T, cache.len, cache.capacity) };
        std::mem::forget(cache);
        v
    }
}

// #[macro_export]
// macro_rules! to_cache {
//     ($v:ident) => {{
//         $v.clear();
//         $v.into_iter().map(|_| unreachable!()).collect_vec().into()
//     }};
// }

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    // function that pushes elements on the cache
    fn use_cache_first(inactivated_cache: CacheMisere<u64>) -> CacheMisere<u64> {
        let mut cache: Vec<&u64> = inactivated_cache.into();
        let data = (0..100).collect_vec();
        assert_eq!(cache.len(), 0);
        assert_eq!(cache.capacity(), 0);

        for i in &data {
            cache.push(i);
        }
        assert_eq!(cache.len(), 100);
        assert_eq!(cache.capacity(), 128);

        cache.into()
    }

    // function that pushes elements on the cache
    // used to detect if the cache retains the capacity
    fn use_cache_second(inactivated_cache: CacheMisere<u64>) -> CacheMisere<u64> {
        let mut cache: Vec<&u64> = inactivated_cache.into();

        let data = (0..100).collect_vec();
        assert_eq!(cache.len(), 0);
        assert_eq!(cache.capacity(), 128); // capacity re used

        for i in &data {
            cache.push(i);
        }
        assert_eq!(cache.len(), 100);
        assert_eq!(cache.capacity(), 128);

        cache.into()
    }

    #[test]
    fn test_cache() {
        let id_bucket: usize = 5;
        let id_hk: usize = 9;
        let is_large: bool = false;
        let cache_value = CachedValue::new(id_bucket, id_hk, is_large);

        assert_eq!(cache_value.get_id_bucket(), id_bucket);
        assert_eq!(cache_value.get_id_hk(), id_hk);
        assert_eq!(cache_value.get_is_large(), is_large);
    }

    #[test]
    fn test_cache_misere() {
        let inactivated_cache = CacheMisere::<u64>::new();
        let inactivated_cache = use_cache_first(inactivated_cache);
        let inactivated_cache = use_cache_second(inactivated_cache);
        drop(inactivated_cache);
    }

    #[test]
    fn test_cache_misere_default() {
        let inactivated_cache: CacheMisere<u64> = Default::default();
        let inactivated_cache = use_cache_first(inactivated_cache);
        let inactivated_cache = use_cache_second(inactivated_cache);
        drop(inactivated_cache);
    }
}
