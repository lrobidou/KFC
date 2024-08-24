use serde::{Deserialize, Serialize};

/// Contains values packed in a single u64 value.
/// Values are:
/// index: the index of a extended hyperkmer in the vector
/// start: the start of the hyperkmer in the extended hyperkmer whose index is `index`
/// end: the end of the hyperkmer in the extended hyperkmer whose index is `index`
/// change_orientation: true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
pub struct HKMetadata {
    data: HKMetadataInner<32, 15, 15>,
}

impl HKMetadata {
    pub fn new(
        index: usize,
        start: usize,
        end: usize,
        is_large: bool,
        // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
        change_orientation: bool,
    ) -> Self {
        Self {
            data: HKMetadataInner::new(index, start, end, is_large, change_orientation),
        }
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

#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, PartialOrd, Ord, Eq)]
pub struct HKMetadataInner<const SIZE_INDEX: usize, const SIZE_START: usize, const SIZE_END: usize>
{
    data: usize,
}

impl<const SIZE_INDEX: usize, const SIZE_START: usize, const SIZE_END: usize>
    HKMetadataInner<SIZE_INDEX, SIZE_START, SIZE_END>
{
    pub fn new(
        index: usize,
        start: usize,
        end: usize,
        is_large: bool,
        // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
        change_orientation: bool,
    ) -> Self {
        debug_assert_eq!(SIZE_INDEX + SIZE_START + SIZE_END + 2, 64);
        debug_assert!(index < 2usize.pow((SIZE_INDEX + 1) as u32));
        debug_assert!(start < 2usize.pow((SIZE_START + 1) as u32));
        debug_assert!(end < 2usize.pow((SIZE_END + 1) as u32));
        let data = index << (SIZE_START + SIZE_END + 2);
        let data = data + ((start << (SIZE_INDEX + SIZE_END + 2)) >> SIZE_INDEX);
        let data = data + ((end << (SIZE_INDEX + SIZE_START + 2)) >> (SIZE_INDEX + SIZE_START));
        let data = data + 2 * (is_large as usize);
        let data = data + change_orientation as usize;

        let me = Self { data };
        assert_eq!(me.get_index(), index);
        assert_eq!(me.get_start(), start);
        assert_eq!(me.get_end(), end);
        assert_eq!(me.get_is_large(), is_large);
        assert_eq!(me.get_change_orientation(), change_orientation);
        me
    }

    pub fn get_index(&self) -> usize {
        // first N bits
        let shift_amount = 64 - SIZE_INDEX;
        self.data >> shift_amount
    }

    pub fn get_start(&self) -> usize {
        let nb_bits_up = SIZE_INDEX;
        let nb_bits_down = SIZE_END + 2;
        (self.data << nb_bits_up) >> (nb_bits_up + nb_bits_down)
    }

    pub fn get_end(&self) -> usize {
        let nb_bits_up = SIZE_INDEX + SIZE_START;
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
