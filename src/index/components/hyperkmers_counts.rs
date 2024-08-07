use super::extended_hyperkmers::ExtendedHyperkmers;
use super::get_subsequence_from_metadata;
use crate::Superkmer;
use crate::{
    check_equal_mashmap,
    superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata},
    Count, Minimizer,
};
use mashmap::MashMap;
use serde::{
    de::{MapAccess, Visitor},
    ser::SerializeMap,
    Deserialize, Deserializer, Serialize, Serializer,
};

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
        // println!("{}", is_large);
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
        // assert_eq!((data >> 1) as bool, is_large);
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
        // last bit of the data
        (self.data >> 1) % 2 == 1
    }

    pub fn get_change_orientation(&self) -> bool {
        // last bit of the data
        self.data % 2 == 1
    }
}

pub struct HKCount {
    data: MashMap<Minimizer, (HKMetadata, HKMetadata, Count)>,
}

impl PartialEq for HKCount {
    fn eq(&self, other: &Self) -> bool {
        check_equal_mashmap(&self.data, &other.data)
    }
}

impl Serialize for HKCount {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.data.len()))?;
        for (k, v) in self.data.iter() {
            map.serialize_entry(k, v)?;
        }
        map.end()
    }
}

struct HKCountVisitor {}

impl HKCountVisitor {
    fn new() -> Self {
        Self {}
    }
}

impl<'de> Visitor<'de> for HKCountVisitor {
    type Value = HKCount;

    // Format a message stating what data this Visitor expects to receive.
    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("hyper kmer counts")
    }

    fn visit_map<M>(self, mut access: M) -> Result<Self::Value, M::Error>
    where
        M: MapAccess<'de>,
    {
        let mut super_kmer_counts_data =
            MashMap::<Minimizer, (HKMetadata, HKMetadata, Count)>::with_capacity(
                access.size_hint().unwrap_or(0),
            );

        while let Some((key, value)) = access.next_entry()? {
            super_kmer_counts_data.insert(key, value);
        }

        Ok(HKCount {
            data: super_kmer_counts_data,
        })
    }
}

impl<'de> Deserialize<'de> for HKCount {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_map(HKCountVisitor::new())
    }
}

impl HKCount {
    pub fn new() -> Self {
        Self {
            data: MashMap::new(),
        }
    }

    pub fn get_data(&self) -> &MashMap<Minimizer, (HKMetadata, HKMetadata, Count)> {
        &self.data
    }

    /// Searches for a match between `extended_hyperkmer_left` and one of the left extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_left_id(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        minimizer: &Minimizer,
        extended_hyperkmer_left: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<(usize, bool)> {
        let canonical_extended_hyperkmer_left = extended_hyperkmer_left.to_canonical();

        for (candidate_left_hk_metadata, _candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let index = candidate_left_hk_metadata.get_index();
            let is_large = candidate_left_hk_metadata.get_is_large();
            let subseq = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_left_hk_metadata,
            );

            if canonical_extended_hyperkmer_left.equal_bitpacked(&subseq) {
                return Some((index, is_large));
            }
        }
        None
    }

    /// Searches for a match between `extended_hyperkmer_right` and one of the right extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_right_id(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        minimizer: &Minimizer,
        extended_hyperkmer_right: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<(usize, bool)> {
        let canonical_extended_hyperkmer_right = extended_hyperkmer_right.to_canonical();
        for (_candidate_left_hk_metadata, candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let index = candidate_right_hk_metadata.get_index();
            let is_large = candidate_right_hk_metadata.get_is_large();

            let subseq = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_right_hk_metadata,
            );

            if canonical_extended_hyperkmer_right.equal_bitpacked(&subseq) {
                return Some((index, is_large));
            }
        }
        None
    }

    /// Add a new entry in the hyperkmer counts.
    /// Minimizer should already be in canonical form.
    pub fn insert_new_entry_in_hyperkmer_count(
        &mut self,
        minimizer: &Minimizer,
        left_metadata: &HKMetadata,
        right_metadata: &HKMetadata,
        count: Count,
    ) {
        assert!(left_metadata.get_start() < left_metadata.get_end());
        assert!(right_metadata.get_start() < right_metadata.get_end());
        self.data
            .insert(*minimizer, (*left_metadata, *right_metadata, count));
    }

    /// Search if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`
    /// If so, increase the count of the occurence and return `true`
    /// Else, return false
    pub fn increase_count_if_exact_match(
        &self,
        minimizer: &Minimizer,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        left_hk: &SubsequenceMetadata<NoBitPacked>,
        right_hk: &SubsequenceMetadata<NoBitPacked>,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
            self.data.get_mut_iter(minimizer)
        {
            let is_exact_match = search_exact_hyperkmer_match(
                hyperkmers,
                large_hyperkmers,
                left_hk,
                right_hk,
                candidate_left_ext_hk_metadata,
                candidate_right_ext_hk_metadata,
            );
            if is_exact_match {
                *count_hk += 1;
                return true;
            }
        }
        false
    }

    pub fn search_for_inclusion(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        superkmer: &Superkmer,
        left_sk: &SubsequenceMetadata<NoBitPacked>,
        right_sk: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let minimizer = superkmer.get_minimizer();
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count_hk) in
            self.data.get_mut_iter(&minimizer)
        {
            // get sequences as they would appear if the current superkmer was canonical
            let is_large_left = candidate_left_ext_hk_metadata.get_is_large();
            let subseq_left = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            );
            let candidate_left_ext_hk: SubsequenceMetadata<BitPacked> = subseq_left
                .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());

            let is_large_right = candidate_right_ext_hk_metadata.get_is_large();
            let subseq_right = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_right_ext_hk_metadata,
            );
            let candidate_right_ext_hk: SubsequenceMetadata<BitPacked> = subseq_right
                .change_orientation_if(candidate_right_ext_hk_metadata.get_change_orientation());

            let match_start_left = candidate_left_ext_hk.starts_with_nobitpacked(left_sk);
            let match_end_left = candidate_left_ext_hk.ends_with_nobitpacked(left_sk);
            let match_left = match_start_left || match_end_left;

            let match_start_right = candidate_right_ext_hk.starts_with_nobitpacked(right_sk);
            let match_end_right = candidate_right_ext_hk.ends_with_nobitpacked(right_sk);
            let match_right = match_start_right || match_end_right;

            if match_left && match_right {
                let (start_left, end_left) = if match_start_left {
                    (0, left_sk.len())
                } else {
                    (
                        candidate_left_ext_hk.len() - left_sk.len(),
                        candidate_left_ext_hk.len(),
                    )
                };

                let (start_right, end_right) = if match_start_right {
                    (0, right_sk.len())
                } else {
                    (
                        candidate_right_ext_hk.len() - right_sk.len(),
                        candidate_right_ext_hk.len(),
                    )
                };

                return Some((
                    HKMetadata::new(
                        candidate_left_ext_hk_metadata.get_index(),
                        start_left,
                        end_left,
                        is_large_left,
                        candidate_left_ext_hk_metadata.get_change_orientation(),
                    ),
                    HKMetadata::new(
                        candidate_right_ext_hk_metadata.get_index(),
                        start_right,
                        end_right,
                        is_large_right,
                        candidate_right_ext_hk_metadata.get_change_orientation(),
                    ),
                ));
            }
        }
        None
    }

    /// Return true if the minimizer is present, false otherwise
    pub fn contains_minimizer(&self, minimizer: &Minimizer) -> bool {
        self.data.contains_key(minimizer)
    }

    /// Search for the maximal inclusion of `left_sk` and `right_sk`
    /// in the extended hyperkmers associated with the minimizers
    /// Returns the metadata associated with this inclusion
    pub fn search_for_maximal_inclusion(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        k: usize,
        m: usize,
        minimizer: &Minimizer,
        left_sk: &SubsequenceMetadata<NoBitPacked>,
        right_sk: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let mut match_size = k - 1;
        let mut match_metadata = None;

        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            // let is_large_left = candidate_left_ext_hk_metadata.get_is_large();
            let subseq_left = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            );
            let candidate_left_ext_hk: SubsequenceMetadata<BitPacked> = subseq_left
                .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());

            // let is_large_right = candidate_right_ext_hk_metadata.get_is_large();
            let subseq_right = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_right_ext_hk_metadata,
            );
            let candidate_right_ext_hk: SubsequenceMetadata<BitPacked> = subseq_right
                .change_orientation_if(candidate_right_ext_hk_metadata.get_change_orientation());

            // extract candidate hyperkmers
            let candidate_left_hyperkmer = &candidate_left_ext_hk.subsequence(
                candidate_left_ext_hk_metadata.get_start(),
                candidate_left_ext_hk_metadata.get_end(),
            );
            let candidate_right_hyperkmer = &candidate_right_ext_hk.subsequence(
                candidate_right_ext_hk_metadata.get_start(),
                candidate_right_ext_hk_metadata.get_end(),
            );

            let len_current_match_left =
                left_sk.common_suffix_length_with_bitpacked(candidate_left_hyperkmer);
            let len_current_match_right =
                right_sk.common_prefix_length_with_bitpacked(candidate_right_hyperkmer);
            let current_match_size = len_current_match_left + len_current_match_right;

            assert!(len_current_match_left >= m - 1);
            assert!(len_current_match_right >= m - 1);

            if current_match_size - 2 * (m - 1) + m > match_size {
                match_size = current_match_size;

                match_metadata = Some((
                    // same suffix => same end, but different start
                    HKMetadata::new(
                        candidate_left_ext_hk_metadata.get_index(),
                        candidate_left_ext_hk_metadata.get_end() - len_current_match_left,
                        candidate_left_ext_hk_metadata.get_end(),
                        candidate_left_ext_hk_metadata.get_is_large(),
                        candidate_left_ext_hk_metadata.get_change_orientation(),
                    ),
                    // same prefix => same start, but different end
                    HKMetadata::new(
                        candidate_right_ext_hk_metadata.get_index(),
                        candidate_right_ext_hk_metadata.get_start(),
                        candidate_right_ext_hk_metadata.get_start() + len_current_match_right,
                        candidate_right_ext_hk_metadata.get_is_large(),
                        candidate_right_ext_hk_metadata.get_change_orientation(),
                    ),
                ));
            }
        }
        match_metadata
    }

    pub fn count_occurence_kmer(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        large_hyperkmers: &Vec<(usize, Vec<u8>)>,
        minimizer: &Minimizer,
        left_context: &SubsequenceMetadata<NoBitPacked>,
        right_context: &SubsequenceMetadata<NoBitPacked>,
        k: usize,
        m: usize,
    ) -> Count {
        let mut total_count = 0;
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count) in
            self.data.get_iter(minimizer)
        {
            // get sequences as they would appear if the current superkmer was canonical
            let candidate_left_hyperkmer = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            )
            .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());
            let candidate_left_hyperkmer = candidate_left_hyperkmer.subsequence(
                candidate_left_ext_hk_metadata.get_start(),
                candidate_left_ext_hk_metadata.get_end(),
            );

            let candidate_right_hyperkmer = get_subsequence_from_metadata(
                hyperkmers,
                large_hyperkmers,
                candidate_right_ext_hk_metadata,
            )
            .change_orientation_if(candidate_right_ext_hk_metadata.get_change_orientation());
            let candidate_right_hyperkmer = candidate_right_hyperkmer.subsequence(
                candidate_right_ext_hk_metadata.get_start(),
                candidate_right_ext_hk_metadata.get_end(),
            );
            let len_current_match_left =
                left_context.common_suffix_length_with_bitpacked(&candidate_left_hyperkmer);
            let len_current_match_right =
                right_context.common_prefix_length_with_bitpacked(&candidate_right_hyperkmer);
            let current_match_size = len_current_match_left + len_current_match_right;

            #[cfg(debug_assertions)]
            {
                let left_string = candidate_left_hyperkmer.to_string();
                let right_string = candidate_right_hyperkmer.to_string();

                debug_assert_eq!(
                    left_string[(left_string.len() - (m - 2))..left_string.len()],
                    right_string[0..(m - 2)]
                );
            }

            if current_match_size - 2 * (m - 1) + m >= k {
                total_count += count;
            }
        }
        total_count
    }
}

pub fn search_exact_hyperkmer_match(
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &Vec<(usize, Vec<u8>)>,
    left_hk: &SubsequenceMetadata<NoBitPacked>,
    right_hk: &SubsequenceMetadata<NoBitPacked>,
    candidate_left_ext_hk_metadata: &HKMetadata,
    candidate_right_ext_hk_metadata: &HKMetadata,
) -> bool {
    // get sequences as they would appear if the current superkmer was canonical
    let candidate_left_hyperkmer =
        get_subsequence_from_metadata(hyperkmers, large_hyperkmers, candidate_left_ext_hk_metadata)
            .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());
    let candidate_left_hyperkmer = candidate_left_hyperkmer.subsequence(
        candidate_left_ext_hk_metadata.get_start(),
        candidate_left_ext_hk_metadata.get_end(),
    );

    let candidate_right_hyperkmer = get_subsequence_from_metadata(
        hyperkmers,
        large_hyperkmers,
        candidate_right_ext_hk_metadata,
    )
    .change_orientation_if(candidate_right_ext_hk_metadata.get_change_orientation());
    let candidate_right_hyperkmer = candidate_right_hyperkmer.subsequence(
        candidate_right_ext_hk_metadata.get_start(),
        candidate_right_ext_hk_metadata.get_end(),
    );

    let match_left = left_hk.equal_bitpacked(&candidate_left_hyperkmer);
    let match_right = right_hk.equal_bitpacked(&candidate_right_hyperkmer);

    match_left && match_right
}
