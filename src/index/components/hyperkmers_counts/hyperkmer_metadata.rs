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
    data: HKMetadataInner<8, 24, 15, 15>,
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
            data: HKMetadataInner::new(bucket_id, index, start, end, is_large, change_orientation),
        }
    }

    pub fn get_bucket_id(&self) -> usize {
        self.data.get_bucket_id()
    }

    pub fn get_index(&self) -> usize {
        self.data.get_index()
    }

    pub fn get_start(&self) -> usize {
        self.data.get_start()
    }

    pub fn get_end(&self) -> usize {
        self.data.get_end()
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
            .field("bucket_id", &self.get_bucket_id())
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
    data: usize,
}

impl<
        const SIZE_BUCKET_ID: usize,
        const SIZE_INDEX: usize,
        const SIZE_START: usize,
        const SIZE_END: usize,
    > HKMetadataInner<SIZE_BUCKET_ID, SIZE_INDEX, SIZE_START, SIZE_END>
{
    pub fn new(
        bucket_id: usize,
        index: usize,
        start: usize,
        end: usize,
        is_large: bool,
        // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
        change_orientation: bool,
    ) -> Self {
        debug_assert_eq!(SIZE_BUCKET_ID + SIZE_INDEX + SIZE_START + SIZE_END + 2, 64);
        assert!(bucket_id < 2usize.pow((SIZE_BUCKET_ID + 1) as u32));
        assert!(index < 2usize.pow((SIZE_INDEX + 1) as u32));
        assert!(start < 2usize.pow((SIZE_START + 1) as u32));
        assert!(end < 2usize.pow((SIZE_END + 1) as u32));
        let data = bucket_id << (SIZE_INDEX + SIZE_START + SIZE_END + 2);
        let data = data + (index << (SIZE_BUCKET_ID + SIZE_START + SIZE_END + 2) >> SIZE_BUCKET_ID);
        let data = data
            + ((start << (SIZE_BUCKET_ID + SIZE_INDEX + SIZE_END + 2))
                >> (SIZE_BUCKET_ID + SIZE_INDEX));
        let data = data
            + ((end << (SIZE_BUCKET_ID + SIZE_INDEX + SIZE_START + 2))
                >> (SIZE_BUCKET_ID + SIZE_INDEX + SIZE_START));
        let data = data + 2 * (is_large as usize);
        let data = data + change_orientation as usize;

        let me = Self { data };

        debug_assert_eq!(me.get_bucket_id(), bucket_id);
        debug_assert_eq!(me.get_index(), index);
        debug_assert_eq!(me.get_start(), start);
        debug_assert_eq!(me.get_end(), end);
        debug_assert_eq!(me.get_is_large(), is_large);
        debug_assert_eq!(me.get_change_orientation(), change_orientation);

        me
    }

    pub fn get_bucket_id(&self) -> usize {
        // first N bits
        let shift_amount = 64 - SIZE_BUCKET_ID;
        self.data >> shift_amount
    }

    pub fn get_index(&self) -> usize {
        let nb_bits_up = SIZE_BUCKET_ID;
        let nb_bits_down = SIZE_START + SIZE_END + 2;
        (self.data << nb_bits_up) >> (nb_bits_up + nb_bits_down)
    }

    pub fn get_start(&self) -> usize {
        let nb_bits_up = SIZE_BUCKET_ID + SIZE_INDEX;
        let nb_bits_down = SIZE_END + 2;
        (self.data << nb_bits_up) >> (nb_bits_up + nb_bits_down)
    }

    pub fn get_end(&self) -> usize {
        let nb_bits_up = SIZE_BUCKET_ID + SIZE_INDEX + SIZE_START;
        let nb_bits_down = 2;
        (self.data << nb_bits_up) >> (nb_bits_up + nb_bits_down)
    }

    pub fn get_is_large(&self) -> bool {
        (self.data >> 1) % 2 == 1
    }

    pub fn get_change_orientation(&self) -> bool {
        // last bit of the data
        self.data % 2 == 1
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
        let hk_metadata_inner = HKMetadataInner::<8, 24, 15, 15>::new(
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
        let hk_metadata_inner = HKMetadataInner::<8, 24, 15, 15>::new(
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
        let bucket_id = 255;
        let index = 16777215;
        let start = 32767;
        let end = 32767;
        let is_large = true;
        let change_orientation = true;

        // creation of the structure
        let hk_metadata_inner = HKMetadataInner::<8, 24, 15, 15>::new(
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
