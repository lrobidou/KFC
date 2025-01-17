use crate::index::components::{HKMetadata, SIZE_BUCKET_ID};

#[derive(Debug, Clone)]
pub struct CachedValue {
    metadata_bucket_index_large: u64,
}

/// Cache of a HKMetadata.
/// The cache is not valid accross reads.
/// Since:
/// - the sequence indicated by the metadata depends on the orientation of the minimizer
/// - we want to use the cache for consecutive minimizers
///
/// => the cache needs to be minimizer-agnostic
/// So the cache represents the hyperkmer as it is in the read
/// => the cache should not be reused across reads
impl CachedValue {
    pub fn from_hk_metadata(metadata: &HKMetadata) -> Self {
        Self {
            metadata_bucket_index_large: metadata.get_bucket_index_large(),
        }
    }

    /// Selects which `HKMetadata` to put in the cache based on `is_sk_canonical_in_the_read`.
    pub fn from_left_and_right(
        is_current_sk_canonical_in_the_read: bool,
        left: &HKMetadata,
        right: &HKMetadata,
    ) -> Self {
        if is_current_sk_canonical_in_the_read {
            Self::from_hk_metadata(right)
        } else {
            Self::from_hk_metadata(left)
        }
    }

    pub fn get_id_bucket(&self) -> usize {
        // first N bits
        let shift_amount = 64 - SIZE_BUCKET_ID;
        (self.metadata_bucket_index_large >> shift_amount)
            .try_into()
            .unwrap()
    }

    pub fn get_id_hk(&self) -> usize {
        let nb_bits_up = SIZE_BUCKET_ID;
        let nb_bits_down = 2;
        ((self.metadata_bucket_index_large << nb_bits_up) >> (nb_bits_up + nb_bits_down))
            .try_into()
            .unwrap()
    }

    pub fn get_is_large(&self) -> bool {
        (self.metadata_bucket_index_large >> 1) % 2 == 1
    }

    // TODO do a full match ?
    pub fn partial_match_metadata(&self, other: &HKMetadata) -> bool {
        self.metadata_bucket_index_large == other.get_bucket_index_large()
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
        let start = 10;
        let end = 56;
        let metadata = HKMetadata::new(id_bucket, id_hk, start, end, is_large, false);
        let cache_value = CachedValue::from_hk_metadata(&metadata);

        assert_eq!(cache_value.get_id_bucket(), id_bucket);
        assert_eq!(cache_value.get_id_hk(), id_hk);
        assert_eq!(cache_value.get_is_large(), is_large);
    }
}
