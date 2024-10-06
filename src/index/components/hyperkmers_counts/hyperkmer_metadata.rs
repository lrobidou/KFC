use serde::{Deserialize, Serialize};

/// Contains values packed in a single u64 value.
///
/// Values are:
/// index: the index of a extended hyperkmer in the vector
/// start: the start of the hyperkmer in the extended hyperkmer whose index is `index`
/// end: the end of the hyperkmer in the extended hyperkmer whose index is `index`
/// change_orientation: true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
#[derive(Clone, Copy, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
pub struct HKMetadata {
    // 2**32 is about 17 times lager than the largest human chomosome
    data: HKMetadataInner<16, 46, 32, 32>,
}

impl HKMetadata {
    pub fn new(
        bucket_id: usize,
        index: usize,
        start: usize,
        end: usize,
        is_large: bool,
        // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
        change_orientation: bool,
    ) -> Self {
        Self {
            data: HKMetadataInner::new(
                bucket_id.try_into().unwrap(),
                index.try_into().unwrap(),
                start.try_into().unwrap(),
                end.try_into().unwrap(),
                is_large,
                change_orientation,
            ),
        }
    }

    pub fn get_bucket_id(&self) -> usize {
        self.data.get_bucket_id().try_into().unwrap()
    }

    pub fn get_index(&self) -> usize {
        self.data.get_index().try_into().unwrap()
    }

    pub fn get_start(&self) -> usize {
        self.data.get_start().try_into().unwrap()
    }

    pub fn get_end(&self) -> usize {
        self.data.get_end().try_into().unwrap()
    }

    pub fn get_is_large(&self) -> bool {
        self.data.get_is_large()
    }

    pub fn get_change_orientation(&self) -> bool {
        self.data.get_change_orientation()
    }
}

impl std::fmt::Debug for HKMetadata {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HKMetadata")
            .field("bucket_id", &self.get_bucket_id())
            .field("index", &self.get_index())
            .field("start", &self.get_start())
            .field("end", &self.get_end())
            .field("is_large", &self.get_is_large())
            .field("change_orientation", &self.get_change_orientation())
            .finish()
    }
}

/// Contains values packed in a single u64 value.
///
/// This is not publicly available, as specifying the size of every field everywhere is cumbersome and error prone.
/// Use `HKMetadata` instead.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
struct HKMetadataInner<
    const SIZE_BUCKET_ID: usize,
    const SIZE_INDEX: usize,
    const SIZE_START: usize,
    const SIZE_END: usize,
> {
    // BUCKET_ID | INDEX | IS_LARGE | CHANGE_ORIENTATION
    bucket_index_end_large: u64,
    // START | END
    start_stop: u64,
}

impl<
        const SIZE_BUCKET_ID: usize,
        const SIZE_INDEX: usize,
        const SIZE_START: usize,
        const SIZE_END: usize,
    > HKMetadataInner<SIZE_BUCKET_ID, SIZE_INDEX, SIZE_START, SIZE_END>
{
    pub fn new(
        bucket_id: u64,
        index: u64,
        start: u64,
        end: u64,
        is_large: bool,
        // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
        change_orientation: bool,
    ) -> Self {
        // check that the sum of size is equal to 64
        debug_assert_eq!(SIZE_BUCKET_ID + SIZE_INDEX + 2, 64);
        debug_assert_eq!(SIZE_START + SIZE_END, 64);
        // check that each field is below the limit
        assert!(bucket_id < 2usize.pow((SIZE_BUCKET_ID + 1) as u32) as u64);
        assert!(index < 2usize.pow((SIZE_INDEX + 1) as u32) as u64);
        assert!(start < 2usize.pow((SIZE_START + 1) as u32) as u64);
        assert!(end < 2usize.pow((SIZE_END + 1) as u32) as u64);

        let bucket_index_end_large = bucket_id << (SIZE_INDEX + 2);
        let bucket_index_end_large =
            bucket_index_end_large + (index << (SIZE_BUCKET_ID + 2) >> SIZE_BUCKET_ID);
        let bucket_index_end_large = bucket_index_end_large + 2 * (is_large as u64);
        let bucket_index_end_large = bucket_index_end_large + change_orientation as u64;

        let start_stop = start << SIZE_END;
        let start_stop = start_stop + ((end << SIZE_START) >> (SIZE_START));

        let me = Self {
            bucket_index_end_large,
            start_stop,
        };

        debug_assert_eq!(me.get_bucket_id(), bucket_id);
        debug_assert_eq!(me.get_index(), index);
        debug_assert_eq!(me.get_start(), start);
        debug_assert_eq!(me.get_end(), end);
        debug_assert_eq!(me.get_is_large(), is_large);
        debug_assert_eq!(me.get_change_orientation(), change_orientation);

        me
    }

    pub fn get_bucket_id(&self) -> u64 {
        // first N bits
        let shift_amount = 64 - SIZE_BUCKET_ID;
        self.bucket_index_end_large >> shift_amount
    }

    pub fn get_index(&self) -> u64 {
        let nb_bits_up = SIZE_BUCKET_ID;
        let nb_bits_down = 2;
        (self.bucket_index_end_large << nb_bits_up) >> (nb_bits_up + nb_bits_down)
    }

    pub fn get_start(&self) -> u64 {
        let nb_bits_down = SIZE_END;
        self.start_stop >> nb_bits_down
    }

    pub fn get_end(&self) -> u64 {
        let nb_bits_up = SIZE_START;
        (self.start_stop << nb_bits_up) >> nb_bits_up
    }

    pub fn get_is_large(&self) -> bool {
        (self.bucket_index_end_large >> 1) % 2 == 1
    }

    pub fn get_change_orientation(&self) -> bool {
        // last bit of the data
        self.bucket_index_end_large % 2 == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyper_kmer_metadata_inner_all_0() {
        // parameters
        let bucket_id = 0;
        let index = 0;
        let start = 0;
        let end = 0;
        let is_large = false;
        let change_orientation = false;

        // creation of the structure
        let hk_metadata_inner = HKMetadataInner::<16, 46, 32, 32>::new(
            bucket_id,
            index,
            start,
            end,
            is_large,
            change_orientation,
        );

        // tests
        assert_eq!(hk_metadata_inner.get_bucket_id(), bucket_id);
        assert_eq!(hk_metadata_inner.get_index(), index);
        assert_eq!(hk_metadata_inner.get_start(), start);
        assert_eq!(hk_metadata_inner.get_end(), end);
        assert_eq!(hk_metadata_inner.get_is_large(), is_large);
        assert_eq!(
            hk_metadata_inner.get_change_orientation(),
            change_orientation
        );
    }
    #[test]
    fn test_hyper_kmer_metadata_inner_all_random_values() {
        // parameters
        let bucket_id = 67;
        let index = 678;
        let start = 8;
        let end = 56;
        let is_large = true;
        let change_orientation = false;

        // creation of the structure
        let hk_metadata_inner = HKMetadataInner::<16, 46, 32, 32>::new(
            bucket_id,
            index,
            start,
            end,
            is_large,
            change_orientation,
        );

        // tests
        assert_eq!(hk_metadata_inner.get_bucket_id(), bucket_id);
        assert_eq!(hk_metadata_inner.get_index(), index);
        assert_eq!(hk_metadata_inner.get_start(), start);
        assert_eq!(hk_metadata_inner.get_end(), end);
        assert_eq!(hk_metadata_inner.get_is_large(), is_large);
        assert_eq!(
            hk_metadata_inner.get_change_orientation(),
            change_orientation
        );
    }

    #[test]
    fn test_hyper_kmer_metadata_inner_all_max() {
        // parameters
        let bucket_id = 65536 - 1;
        let index = 70368744177664 - 1;
        let start = 4294967296 - 1;
        let end = 4294967296 - 1;
        let is_large = true;
        let change_orientation = true;

        // creation of the structure
        let hk_metadata_inner = HKMetadataInner::<16, 46, 32, 32>::new(
            bucket_id,
            index,
            start,
            end,
            is_large,
            change_orientation,
        );

        // tests
        assert_eq!(hk_metadata_inner.get_bucket_id(), bucket_id);
        assert_eq!(hk_metadata_inner.get_index(), index);
        assert_eq!(hk_metadata_inner.get_start(), start);
        assert_eq!(hk_metadata_inner.get_end(), end);
        assert_eq!(hk_metadata_inner.get_is_large(), is_large);
        assert_eq!(
            hk_metadata_inner.get_change_orientation(),
            change_orientation
        );
    }
}
