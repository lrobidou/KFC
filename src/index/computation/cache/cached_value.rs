use crate::index::components::HKMetadata;

// todo: use all HKMetadata fields ?
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

    /// Selects which `HKMetadata` to put in the cache based on `is_current_sk_canonical_in_the_read`.
    pub fn from_left_and_right(
        is_current_sk_canonical_in_the_read: bool,
        left: &HKMetadata,
        right: &HKMetadata,
    ) -> Self {
        if is_current_sk_canonical_in_the_read {
            Self::new(
                right.get_bucket_id(),
                right.get_index(),
                right.get_is_large(),
            )
        } else {
            Self::new(left.get_bucket_id(), left.get_index(), left.get_is_large())
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

#[cfg(test)]
mod tests {

    use super::*;

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
}
