mod hyperkmer_metadata;

use itertools::Itertools;
use mashmap::MashMap;
use serde::{
    de::{MapAccess, Visitor},
    ser::SerializeMap,
    Deserialize, Deserializer, Serialize, Serializer,
};
use std::collections::HashSet;
use std::sync::atomic::Ordering::SeqCst;

use super::hyperkmer_parts::HyperkmerParts;
use crate::{
    index::computation::{CacheMisere, CachedValue},
    subsequence::{BitPacked, NoBitPacked, Subsequence},
    AtomicCount, Count, Minimizer, Superkmer,
};

pub use hyperkmer_metadata::HKMetadata;

// Notation: in this module, variables starting with `c_` are `candidate`s variables,
// e.g. variables computed by iterating over the `HKCount` table.

pub struct HKCount {
    data: MashMap<Minimizer, (HKMetadata, HKMetadata, AtomicCount)>,
}

impl PartialEq for HKCount {
    fn eq(&self, other: &Self) -> bool {
        check_equal_mashmap_atomic(&self.data, &other.data)
    }
}

// TODO implement this in mashmap ?
/// Checks if two `MashMap` are equal
fn check_equal_mashmap_atomic<K: Ord>(
    map0: &MashMap<K, (HKMetadata, HKMetadata, AtomicCount)>,
    map1: &MashMap<K, (HKMetadata, HKMetadata, AtomicCount)>,
) -> bool {
    // OPTIMIZE this is a naive implementation
    let mut v0 = map0
        .iter()
        .map(|(k, v)| (k, v.0, v.1, v.2.load(SeqCst)))
        .collect_vec();
    let mut v1 = map1
        .iter()
        .map(|(k, v)| (k, v.0, v.1, v.2.load(SeqCst)))
        .collect_vec();
    v0.sort();
    v1.sort();
    v0 == v1
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

pub enum ExactMatchOrInclusion {
    ExactMatch(HKMetadata, HKMetadata),
    Inclusion(HKMetadata, HKMetadata),
    NotFound,
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
        let mut hyper_kmer_counts_data =
            MashMap::<Minimizer, (HKMetadata, HKMetadata, AtomicCount)>::with_capacity(
                access.size_hint().unwrap_or(0),
            );

        while let Some((key, value)) = access.next_entry()? {
            hyper_kmer_counts_data.insert(key, value);
        }

        Ok(HKCount {
            data: hyper_kmer_counts_data,
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

    pub fn get_data(&self) -> &MashMap<Minimizer, (HKMetadata, HKMetadata, AtomicCount)> {
        &self.data
    }

    pub fn minimizer_set(&self) -> HashSet<u64> {
        let mut set = HashSet::with_capacity(self.data.len());

        for (minimizer, _v) in self.data.iter() {
            set.insert(*minimizer);
        }

        set
    }

    /// Searches for a match between `extended_hyperkmer_left` and one of the left extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_left_id(
        &self,
        hyperkmers: &HyperkmerParts,
        minimizer: &Minimizer,
        extended_hyperkmer_left: &Subsequence<NoBitPacked>,
    ) -> Option<(usize, usize, bool)> {
        let canonical_extended_hyperkmer_left = extended_hyperkmer_left.to_canonical();

        for (candidate_left_hk_metadata, _candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let bucket_id = candidate_left_hk_metadata.get_bucket_id();
            let index = candidate_left_hk_metadata.get_index();
            let is_large = candidate_left_hk_metadata.get_is_large();
            let subseq = hyperkmers.get_subsequence_from_metadata(candidate_left_hk_metadata);

            if canonical_extended_hyperkmer_left.equal_bitpacked(&subseq) {
                return Some((bucket_id, index, is_large));
            }
        }
        None
    }

    /// Searches for a match between `extended_hyperkmer_right` and one of the right extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_right_id(
        &self,
        hyperkmers: &HyperkmerParts,
        minimizer: &Minimizer,
        extended_hyperkmer_right: &Subsequence<NoBitPacked>,
    ) -> Option<(usize, usize, bool)> {
        let canonical_extended_hyperkmer_right = extended_hyperkmer_right.to_canonical();
        for (_candidate_left_hk_metadata, candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let bucket_id = candidate_right_hk_metadata.get_bucket_id();
            let index = candidate_right_hk_metadata.get_index();
            let is_large = candidate_right_hk_metadata.get_is_large();

            let subseq = hyperkmers.get_subsequence_from_metadata(candidate_right_hk_metadata);

            if canonical_extended_hyperkmer_right.equal_bitpacked(&subseq) {
                return Some((bucket_id, index, is_large));
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
        self.data.insert(
            *minimizer,
            (*left_metadata, *right_metadata, AtomicCount::new(count)),
        );
    }

    // /// Search if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`
    // /// If so, increase the count of the occurence and return `true`
    // /// Else, return false
    // pub fn increase_count_if_exact_match(
    //     &self,
    //     minimizer: &Minimizer,
    //     hyperkmers: &HyperkmerParts,
    //     left_hk: &Subsequence<NoBitPacked>,
    //     right_hk: &Subsequence<NoBitPacked>,
    // ) -> bool {
    //     for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
    //         //DEBUG why not self mut ?
    //         self.data.get_iter(minimizer)
    //     {
    //         let is_exact_match = search_exact_hyperkmer_match(
    //             hyperkmers,
    //             left_hk,
    //             right_hk,
    //             candidate_left_ext_hk_metadata,
    //             candidate_right_ext_hk_metadata,
    //         );
    //         if is_exact_match {
    //             count_hk.fetch_add(1, SeqCst);
    //             return true;
    //         }
    //     }
    //     false
    // }

    fn increase_count_if_exact_match_in_cache_only<'a>(
        &'a self,
        minimizer: &'a Minimizer,
        is_minimizer_canonical_in_the_read: bool,
        hyperkmers: &HyperkmerParts,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
        cache_vec: &mut Vec<&'a (HKMetadata, HKMetadata, AtomicCount)>,
        cached_value: &CachedValue,
    ) -> Option<(HKMetadata, HKMetadata)> {
        if is_minimizer_canonical_in_the_read {
            // the cache is for the left
            for entry in self.data.get_iter(minimizer) {
                let (c_left_ext_hk_metadata, c_right_ext_hk_metadata, count_hk) = entry;
                if c_left_ext_hk_metadata.get_bucket_id() == cached_value.get_id_bucket()
                    && c_left_ext_hk_metadata.get_is_large() == cached_value.get_is_large()
                    && c_left_ext_hk_metadata.get_index() == cached_value.get_id_hk()
                    && search_exact_hyperkmer_match_right_first(
                        hyperkmers,
                        left_sk,
                        right_sk,
                        c_left_ext_hk_metadata,
                        c_right_ext_hk_metadata,
                    )
                {
                    count_hk.fetch_add(1, SeqCst);
                    return Some((entry.0, entry.1));
                }
                // not an exact match: if no exact match is going to be found, we will iterate over this minimizer again
                // to prevent iterating over the whole map again, I'll store the entries in a cache and iterate over that cache again
                cache_vec.push(entry);
            }
            None
        } else {
            // the cache is for the right
            for entry in self.data.get_iter(minimizer) {
                let (c_left_ext_hk_metadata, c_right_ext_hk_metadata, count_hk) = entry;
                if c_right_ext_hk_metadata.get_bucket_id() == cached_value.get_id_bucket()
                    && c_right_ext_hk_metadata.get_is_large() == cached_value.get_is_large()
                    && c_right_ext_hk_metadata.get_index() == cached_value.get_id_hk()
                    && search_exact_hyperkmer_match_left_first(
                        hyperkmers,
                        left_sk,
                        right_sk,
                        c_left_ext_hk_metadata,
                        c_right_ext_hk_metadata,
                    )
                {
                    count_hk.fetch_add(1, SeqCst);
                    return Some((entry.0, entry.1));
                }
                // not an exact match: if no exact match is going to be found, we will iterate over this minimizer again
                // to prevent iterating over the whole map again, I'll store the entries in a cache and iterate over that cache again
                cache_vec.push(entry);
            }
            None
        }
    }

    fn search_for_exact_match_and_fill_cache<'a>(
        &'a self,
        minimizer: &'a Minimizer,
        is_minimizer_canonical_in_the_read: bool,
        hyperkmers: &HyperkmerParts,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
        cache_vec: &mut Vec<&'a (HKMetadata, HKMetadata, AtomicCount)>,
        cached_value: &Option<CachedValue>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        match cached_value {
            Some(cached_value) => {
                // strategy: use the cached value to do a first pass
                // if it is not there, re do a pass searching for exact match on all entry

                // first pass in question
                let exact_match = self.increase_count_if_exact_match_in_cache_only(
                    minimizer,
                    is_minimizer_canonical_in_the_read,
                    hyperkmers,
                    left_sk,
                    right_sk,
                    cache_vec,
                    cached_value,
                );
                if exact_match.is_some() {
                    return exact_match;
                }
                // not found by looking at entry equal to the cache only => we need another pass
                for (c_left_ext_hk_metadata, c_right_ext_hk_metadata, count_hk) in cache_vec.iter()
                {
                    let is_exact_match = search_exact_hyperkmer_match(
                        hyperkmers,
                        left_sk,
                        right_sk,
                        c_left_ext_hk_metadata,
                        c_right_ext_hk_metadata,
                    );
                    if is_exact_match {
                        count_hk.fetch_add(1, SeqCst);
                        return Some((*c_left_ext_hk_metadata, *c_right_ext_hk_metadata));
                    }
                }
            }
            None => {
                for entry in self.data.get_iter(minimizer) {
                    let (c_left_ext_hk_metadata, c_right_ext_hk_metadata, count_hk) = entry;

                    let is_exact_match = search_exact_hyperkmer_match(
                        hyperkmers,
                        left_sk,
                        right_sk,
                        c_left_ext_hk_metadata,
                        c_right_ext_hk_metadata,
                    );
                    if is_exact_match {
                        count_hk.fetch_add(1, SeqCst);
                        return Some((*c_left_ext_hk_metadata, *c_right_ext_hk_metadata));
                    }
                    // not an exact match: if no exact match is going to be found, we will iterate over this minimizer again
                    // to prevent iterating over the whole map again, I'll store the entries in a cache and iterate over that cache again
                    cache_vec.push(entry);
                }
            }
        }

        #[cfg(debug_assertions)]
        {
            let expected_cache = self.data.get_iter(minimizer).collect_vec();
            assert_eq!(expected_cache.len(), cache_vec.len());
            for (x, y) in cache_vec.iter().zip(expected_cache) {
                assert_eq!(x.0, y.0);
                assert_eq!(x.1, y.1);
                assert_eq!(x.2.load(SeqCst), y.2.load(SeqCst));
            }
        }

        None
    }

    /// Search if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`
    /// If so, increase the count of the occurence and return `true`
    /// Else, return false
    pub fn increase_count_if_exact_match_else_search_for_inclusion(
        &self,
        minimizer: &Minimizer,
        is_minimizer_canonical_in_the_read: bool,
        hyperkmers: &HyperkmerParts,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
        cache: CacheMisere<(HKMetadata, HKMetadata, AtomicCount)>,
        cached_value: &Option<CachedValue>,
    ) -> (
        ExactMatchOrInclusion,
        CacheMisere<(HKMetadata, HKMetadata, AtomicCount)>,
    ) {
        let mut cache_vec: Vec<&(HKMetadata, HKMetadata, AtomicCount)> = cache.into();

        let included = self.search_for_exact_match_and_fill_cache(
            minimizer,
            is_minimizer_canonical_in_the_read,
            hyperkmers,
            left_sk,
            right_sk,
            &mut cache_vec,
            cached_value,
        );

        if let Some(matching_entry) = included {
            return (
                ExactMatchOrInclusion::ExactMatch(matching_entry.0, matching_entry.1),
                cache_vec.into(),
            );
        }

        // no exact match found => we search for inclusion in the cache
        for (c_left_hk_metadata, c_right_hk_metadata, _count_hk) in cache_vec.iter() {
            // extract relevant subsequences frome whole context
            let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                hyperkmers,
                c_left_hk_metadata,
                c_right_hk_metadata,
            );

            let match_start_left = c_left_hk.starts_with_nobitpacked(left_sk);
            let match_end_left = c_left_hk.ends_with_nobitpacked(left_sk);
            let match_left = match_start_left || match_end_left;

            let match_start_right = c_right_hk.starts_with_nobitpacked(right_sk);
            let match_end_right = c_right_hk.ends_with_nobitpacked(right_sk);
            let match_right = match_start_right || match_end_right;

            if match_left && match_right {
                let (start_left, end_left) = if match_start_left {
                    (0, left_sk.len())
                } else {
                    (c_left_hk.len() - left_sk.len(), c_left_hk.len())
                };

                let (start_right, end_right) = if match_start_right {
                    (0, right_sk.len())
                } else {
                    (c_right_hk.len() - right_sk.len(), c_right_hk.len())
                };

                let inclusion = ExactMatchOrInclusion::Inclusion(
                    HKMetadata::new(
                        c_left_hk_metadata.get_bucket_id(),
                        c_left_hk_metadata.get_index(),
                        start_left,
                        end_left,
                        c_left_hk_metadata.get_is_large(),
                        c_left_hk_metadata.get_change_orientation(),
                    ),
                    HKMetadata::new(
                        c_right_hk_metadata.get_bucket_id(),
                        c_right_hk_metadata.get_index(),
                        start_right,
                        end_right,
                        c_right_hk_metadata.get_is_large(),
                        c_right_hk_metadata.get_change_orientation(),
                    ),
                );
                return (inclusion, cache_vec.into());
            }
        }
        (ExactMatchOrInclusion::NotFound, cache_vec.into())
    }

    /// Searches if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`.
    #[cfg(debug_assertions)]
    pub fn search_exact_match(
        &mut self,
        minimizer: &Minimizer,
        hyperkmers: &HyperkmerParts,
        left_hk: &Subsequence<NoBitPacked>,
        right_hk: &Subsequence<NoBitPacked>,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let is_exact_match = search_exact_hyperkmer_match(
                hyperkmers,
                left_hk,
                right_hk,
                candidate_left_ext_hk_metadata,
                candidate_right_ext_hk_metadata,
            );
            if is_exact_match {
                return true;
            }
        }
        false
    }

    pub fn search_for_inclusion(
        &self,
        hyperkmers: &HyperkmerParts,
        superkmer: &Superkmer,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let minimizer = superkmer.get_minimizer();
        for (c_left_hk_metadata, c_right_hk_metadata, _count_hk) in self.data.get_iter(&minimizer) {
            // extract relevant subsequences frome whole context
            let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                hyperkmers,
                c_left_hk_metadata,
                c_right_hk_metadata,
            );

            let match_start_left = c_left_hk.starts_with_nobitpacked(left_sk);
            let match_end_left = c_left_hk.ends_with_nobitpacked(left_sk);
            let match_left = match_start_left || match_end_left;

            let match_start_right = c_right_hk.starts_with_nobitpacked(right_sk);
            let match_end_right = c_right_hk.ends_with_nobitpacked(right_sk);
            let match_right = match_start_right || match_end_right;

            if match_left && match_right {
                let (start_left, end_left) = if match_start_left {
                    (0, left_sk.len())
                } else {
                    (c_left_hk.len() - left_sk.len(), c_left_hk.len())
                };

                let (start_right, end_right) = if match_start_right {
                    (0, right_sk.len())
                } else {
                    (c_right_hk.len() - right_sk.len(), c_right_hk.len())
                };

                return Some((
                    HKMetadata::new(
                        c_left_hk_metadata.get_bucket_id(),
                        c_left_hk_metadata.get_index(),
                        start_left,
                        end_left,
                        c_left_hk_metadata.get_is_large(),
                        c_left_hk_metadata.get_change_orientation(),
                    ),
                    HKMetadata::new(
                        c_right_hk_metadata.get_bucket_id(),
                        c_right_hk_metadata.get_index(),
                        start_right,
                        end_right,
                        c_right_hk_metadata.get_is_large(),
                        c_right_hk_metadata.get_change_orientation(),
                    ),
                ));
            }
        }
        None
    }

    /// Return true if the minimizer is present, false otherwise
    #[inline]
    pub fn contains_minimizer(&self, minimizer: &Minimizer) -> bool {
        self.data.contains_key(minimizer)
    }

    /// Search for the maximal inclusion of `left_sk` and `right_sk`
    /// in the extended hyperkmers associated with the minimizers
    /// Returns the metadata associated with this inclusion
    pub fn search_for_maximal_inclusion(
        &self,
        hyperkmers: &HyperkmerParts,
        k: usize,
        m: usize,
        minimizer: &Minimizer,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let mut match_size = k - 1;
        let mut match_metadata = None;
        // maximal posssible inclusion: if we reach this, we can exit early
        // caution: this includes the (m - 1) duplicationbases in the minimizer (twice)
        let max_inclusion = left_sk.len() + right_sk.len();

        for (c_left_hk_metadata, c_right_hk_metadata, _count) in self.data.get_iter(minimizer) {
            // get sequences as they would appear if the current superkmer was canonical
            // extract relevant subsequences frome whole context
            let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                hyperkmers,
                c_left_hk_metadata,
                c_right_hk_metadata,
            );

            let len_current_match_left = left_sk.common_suffix_length_with_bitpacked(&c_left_hk);
            let len_current_match_right = right_sk.common_prefix_length_with_bitpacked(&c_right_hk);
            let current_match_size = len_current_match_left + len_current_match_right;

            debug_assert!(len_current_match_left >= m - 1);
            debug_assert!(len_current_match_right >= m - 1);

            if current_match_size - 2 * (m - 1) + m > match_size {
                match_size = current_match_size;

                match_metadata = Some((
                    // same suffix => same end, but different start
                    HKMetadata::new(
                        c_left_hk_metadata.get_bucket_id(),
                        c_left_hk_metadata.get_index(),
                        c_left_hk_metadata.get_end() - len_current_match_left,
                        c_left_hk_metadata.get_end(),
                        c_left_hk_metadata.get_is_large(),
                        c_left_hk_metadata.get_change_orientation(),
                    ),
                    // same prefix => same start, but different end
                    HKMetadata::new(
                        c_right_hk_metadata.get_bucket_id(),
                        c_right_hk_metadata.get_index(),
                        c_right_hk_metadata.get_start(),
                        c_right_hk_metadata.get_start() + len_current_match_right,
                        c_right_hk_metadata.get_is_large(),
                        c_right_hk_metadata.get_change_orientation(),
                    ),
                ));

                debug_assert!(current_match_size <= max_inclusion);

                // nothing could go beyond
                if current_match_size == max_inclusion {
                    return match_metadata;
                }
            }
        }
        match_metadata
    }

    pub fn count_occurence_kmer(
        &self,
        hyperkmers: &HyperkmerParts,
        minimizer: &Minimizer,
        left_context: &Subsequence<NoBitPacked>,
        right_context: &Subsequence<NoBitPacked>,
        k: usize,
        m: usize,
    ) -> Count {
        let mut total_count = 0;
        for (c_left_hk_metadata, c_right_hk_metadata, count) in self.data.get_iter(minimizer) {
            // extract relevant subsequences frome whole context
            let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                hyperkmers,
                c_left_hk_metadata,
                c_right_hk_metadata,
            );

            let len_current_match_left =
                left_context.common_suffix_length_with_bitpacked(&c_left_hk);
            let len_current_match_right =
                right_context.common_prefix_length_with_bitpacked(&c_right_hk);
            let current_match_size = len_current_match_left + len_current_match_right;

            #[cfg(debug_assertions)]
            {
                let left_string = c_left_hk.to_string();
                let right_string = c_right_hk.to_string();

                debug_assert_eq!(
                    left_string[(left_string.len() - (m - 2))..left_string.len()],
                    right_string[0..(m - 2)]
                );
            }

            if current_match_size - 2 * (m - 1) + m >= k {
                total_count += count.load(SeqCst);
            }
        }
        total_count
    }

    /// Finds a context that is an exact match on the right and an inclusion on the left. If this context exists, increments the count.
    pub fn increase_count_of_sk_if_found_else_insert_it(
        &mut self,
        k: usize,
        hyperkmers: &HyperkmerParts,
        minimizer: &Minimizer,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
    ) {
        let mut match_case = MatchCases::default();

        // let mut found_match_right = false;
        for (c_left_hk_metadata, c_right_hk_metadata, count) in self.data.get_mut_iter(minimizer) {
            // extract relevant subsequences frome whole context
            let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                hyperkmers,
                c_left_hk_metadata,
                c_right_hk_metadata,
            );

            let left_suffix_size = left_sk.common_suffix_length_with_bitpacked(&c_left_hk);
            let right_prefix_size = right_sk.common_prefix_length_with_bitpacked(&c_right_hk);

            let left_match = if left_suffix_size == left_sk.len() {
                if left_sk.len() == c_left_hk.len() {
                    MatchCase::Exact(*c_left_hk_metadata)
                } else {
                    let new_left_metadata = HKMetadata::new(
                        c_left_hk_metadata.get_bucket_id(),
                        c_left_hk_metadata.get_index(),
                        c_left_hk_metadata.get_end() - left_suffix_size,
                        c_left_hk_metadata.get_end(),
                        c_left_hk_metadata.get_is_large(),
                        c_left_hk_metadata.get_change_orientation(),
                    );
                    MatchCase::Inclusion(new_left_metadata)
                }
            } else {
                // left_prefix_size < left_sk.len()
                MatchCase::Nothing
            };
            let right_match = if right_prefix_size == right_sk.len() {
                if right_sk.len() == c_right_hk.len() {
                    MatchCase::Exact(*c_right_hk_metadata)
                } else {
                    let new_right_metadata = HKMetadata::new(
                        c_right_hk_metadata.get_bucket_id(),
                        c_right_hk_metadata.get_index(),
                        c_right_hk_metadata.get_start(),
                        c_right_hk_metadata.get_start() + right_prefix_size,
                        c_right_hk_metadata.get_is_large(),
                        c_right_hk_metadata.get_change_orientation(),
                    );
                    MatchCase::Inclusion(new_right_metadata)
                }
            } else {
                // right_prefix_size < right_sk.len()
                MatchCase::Nothing
            };

            match_case = std::cmp::min(match_case, MatchCases::new(left_match, right_match));

            // if the cost of the case is low enough, we stop
            // higher threshold => we accept worst solutions (potentially leading to more memory consumption)
            // but this let us stop the search quicker than when using low thresholds => indexation is faster
            let cost = match_case.cost();
            if cost == 0 {
                count.fetch_add(1, SeqCst);
                break;
            } else if match_case.cost() <= 8 {
                break;
            }
        }

        // we identified which match we have
        // let's handle the corresponding case
        match match_case.data() {
            (MatchCase::Exact(_), MatchCase::Exact(_)) => {
                // we do nothing, as the count was already increased
            }
            (MatchCase::Exact(left_hk_metadata), MatchCase::Inclusion(right_hk_metadata)) => {
                // we have to add a `HKMetadata` on the right with a count of 1
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    left_hk_metadata,
                    right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Exact(left_hk_metadata), MatchCase::Nothing) => {
                let right_hk_metadata = insert_new_right_hyperkmer_and_compute_associated_metadata(
                    right_sk, hyperkmers, k,
                );
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    left_hk_metadata,
                    &right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Inclusion(left_hk_metadata), MatchCase::Exact(right_hk_metadata)) => {
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    left_hk_metadata,
                    right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Inclusion(left_hk_metadata), MatchCase::Inclusion(right_hk_metadata)) => {
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    left_hk_metadata,
                    right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Inclusion(left_hk_metadata), MatchCase::Nothing) => {
                let right_hk_metadata = insert_new_right_hyperkmer_and_compute_associated_metadata(
                    right_sk, hyperkmers, k,
                );
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    left_hk_metadata,
                    &right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Nothing, MatchCase::Exact(right_hk_metadata)) => {
                let left_hk_metadata = insert_new_left_hyperkmer_and_compute_associated_metadata(
                    left_sk, hyperkmers, k,
                );
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    &left_hk_metadata,
                    right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Nothing, MatchCase::Inclusion(right_hk_metadata)) => {
                let left_hk_metadata = insert_new_left_hyperkmer_and_compute_associated_metadata(
                    left_sk, hyperkmers, k,
                );
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    &left_hk_metadata,
                    right_hk_metadata,
                    1,
                );
            }
            (MatchCase::Nothing, MatchCase::Nothing) => {
                let left_hk_metadata = insert_new_left_hyperkmer_and_compute_associated_metadata(
                    left_sk, hyperkmers, k,
                );
                let right_hk_metadata = insert_new_right_hyperkmer_and_compute_associated_metadata(
                    right_sk, hyperkmers, k,
                );
                self.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    &left_hk_metadata,
                    &right_hk_metadata,
                    1,
                );
            }
        };
    }
}

/// Create a new (not large) hyperkmer and returns the associated `HKMetadata`.
fn insert_new_right_hyperkmer_and_compute_associated_metadata(
    sequence: &Subsequence<NoBitPacked>,
    hyperkmers: &HyperkmerParts,
    k: usize,
) -> HKMetadata {
    let nb_base_missing = (k - 1) - sequence.len();
    let hyperkmers = hyperkmers.get_typical_parts();
    let metadata = if sequence.is_canonical() {
        assert!(sequence.is_canonical());
        let mut start = sequence.as_vec();
        start.extend_from_slice(&vec![b'A'; nb_base_missing]);
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);

        #[cfg(debug_assertions)]
        {
            // Adding A at the end of a canonical read keeps it canonical
            if !full_sequence.is_canonical() {
                // oops, we juste changed the canonicity of this sequence...
                // ... or did we ?
                // maybe it was a palindrome all along !
                if !full_sequence.is_equal_to_its_revcomp() {
                    panic!("Adding A at the end of a canonical read made it canonical");
                }
            }
        }

        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(&full_sequence);
        // create the associated metadata

        HKMetadata::new(id_bucket, id_hk, 0, sequence.len(), false, false)
    } else {
        let mut start = sequence.as_vec();
        start.extend_from_slice(&vec![b'T'; nb_base_missing]);
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        // TODO prove it

        #[cfg(debug_assertions)]
        {
            // Adding T at the end of a non canonical read does not make it canonical
            if full_sequence.is_canonical() {
                // oops, we juste changed the canonicity of this sequence...
                // ... or did we ?
                // maybe it was a palindrome all along !
                if !full_sequence.is_equal_to_its_revcomp() {
                    panic!();
                }
            }
        }
        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(&full_sequence);
        // create the associated metadata
        HKMetadata::new(id_bucket, id_hk, 0, sequence.len(), false, true)
    };
    metadata
}

/// Create a new (not large) hyperkmer and returns the associated `HKMetadata`.
/// Since this is the left hyperkmer, if it is less the k-1 bases, we need to add some base at its beginning
fn insert_new_left_hyperkmer_and_compute_associated_metadata(
    sequence: &Subsequence<NoBitPacked>,
    hyperkmers: &HyperkmerParts,
    k: usize,
) -> HKMetadata {
    let nb_base_missing = (k - 1) - sequence.len();
    let hyperkmers = hyperkmers.get_typical_parts();
    if sequence.is_canonical() {
        assert!(sequence.is_canonical());
        let mut start = vec![b'A'; nb_base_missing];
        start.extend_from_slice(&sequence.as_vec());
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        assert!(full_sequence.is_canonical()); // Adding A at the beginning of a canonical read keep it canonical

        #[cfg(debug_assertions)]
        {
            // Adding T at the beginning of a non canonical read does not make it canonical
            if !full_sequence.is_canonical() {
                // oops, we juste changed the canonicity of this sequence...
                // ... or did we ?
                // maybe it was a palindrome all along !
                if !full_sequence.is_equal_to_its_revcomp() {
                    panic!();
                }
            }
        }

        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(&full_sequence);
        // create the associated metadata

        HKMetadata::new(
            id_bucket,
            id_hk,
            (k - 1) - sequence.len(),
            k - 1,
            false,
            false,
        )
    } else {
        let mut start = vec![b'T'; nb_base_missing];
        start.extend_from_slice(&sequence.as_vec());
        let full_sequence = Subsequence::new(&start, 0, k - 1, true); // DEBUG
        #[cfg(debug_assertions)]
        {
            // Adding T at the beginning of a non canonical read does not make it canonical
            if full_sequence.is_canonical() {
                // oops, we juste changed the canonicity of this sequence...
                // ... or did we ?
                // maybe it was a palindrome all along !
                if !full_sequence.is_equal_to_its_revcomp() {
                    panic!();
                }
            }
        }
        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(&full_sequence);
        // create the associated metadata
        HKMetadata::new(
            id_bucket,
            id_hk,
            (k - 1) - sequence.len(),
            k - 1,
            false,
            true,
        )
    }
}

#[derive(Eq, PartialEq, Debug)]
enum MatchCase {
    Exact(HKMetadata),
    Inclusion(HKMetadata),
    Nothing,
}

#[derive(Eq, PartialEq, Debug)]
struct MatchCases {
    left: MatchCase,
    right: MatchCase,
}

// "cost" tab:
//
//       |              |               left
//----------------------|----------------------------------------------
//       |              | exact match        inclusion       nothing
//-------|--------------|----------------------------------------------
// right |  exact match |     0                  1              4
//       |              |
//       |  inclusion   |     2                  3              6
//       |              |
//       |  nothing     |     5                  7              8
impl MatchCases {
    fn new(left: MatchCase, right: MatchCase) -> Self {
        Self { left, right }
    }

    fn data(&self) -> (&MatchCase, &MatchCase) {
        (&self.left, &self.right)
    }

    fn cost(&self) -> u8 {
        match (&self.left, &self.right) {
            (MatchCase::Exact(_), MatchCase::Exact(_)) => 0,
            (MatchCase::Exact(_), MatchCase::Inclusion(_)) => 2,
            (MatchCase::Exact(_), MatchCase::Nothing) => 5,
            (MatchCase::Inclusion(_), MatchCase::Exact(_)) => 1,
            (MatchCase::Inclusion(_), MatchCase::Inclusion(_)) => 3,
            (MatchCase::Inclusion(_), MatchCase::Nothing) => 7,
            (MatchCase::Nothing, MatchCase::Exact(_)) => 4,
            (MatchCase::Nothing, MatchCase::Inclusion(_)) => 6,
            (MatchCase::Nothing, MatchCase::Nothing) => 8,
        }
    }
}

impl Default for MatchCases {
    fn default() -> Self {
        Self {
            left: MatchCase::Nothing,
            right: MatchCase::Nothing,
        }
    }
}
impl PartialOrd for MatchCases {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MatchCases {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cost().cmp(&other.cost())
    }
}

pub fn extract_subsequence_from_hk_metadata<'a>(
    hyperkmers: &'a HyperkmerParts,
    hk_metadata: &HKMetadata,
) -> Subsequence<BitPacked<'a>> {
    hyperkmers
        .get_subsequence_from_metadata(hk_metadata)
        .change_orientation_if(hk_metadata.get_change_orientation())
        .subsequence(hk_metadata.get_start(), hk_metadata.get_end())
}

pub fn extract_left_and_right_subsequences<'a>(
    hyperkmers: &'a HyperkmerParts,
    left_hk_metadata: &HKMetadata,
    right_hk_metadata: &HKMetadata,
) -> (Subsequence<BitPacked<'a>>, Subsequence<BitPacked<'a>>) {
    (
        extract_subsequence_from_hk_metadata(hyperkmers, left_hk_metadata),
        extract_subsequence_from_hk_metadata(hyperkmers, right_hk_metadata),
    )
}

/// Checks that `left_hk` and `right_hk` are equal to the subsequences stored in hyperkmers
/// indicated by `left_ext_hk_metadata` and `right_ext_hk_metadata`.
pub fn search_exact_hyperkmer_match(
    hyperkmers: &HyperkmerParts,
    left_hk: &Subsequence<NoBitPacked>,
    right_hk: &Subsequence<NoBitPacked>,
    left_ext_hk_metadata: &HKMetadata,
    right_ext_hk_metadata: &HKMetadata,
) -> bool {
    // Maybe the best `search_exact_hyperkmer_match` was the `search_exact_hyperkmer_match_left_first` we met along the way...
    // ...
    // (... because `search_exact_hyperkmer_match_left_first` dereferences only half of the parts
    // (namley: the left one...)
    // if the left does not match,
    // yielding a small performance boost)
    search_exact_hyperkmer_match_left_first(
        hyperkmers,
        left_hk,
        right_hk,
        left_ext_hk_metadata,
        right_ext_hk_metadata,
    )

    // old code (using `search_exact_hyperkmer_match_left_first` gives a ~5% boost on the whole program)
    // (this was measured once on my laptop, so not the best benchmark ever...)
    // OPTIMIZE: bench the two approaches correctly
    // // extract relevant subsequences frome whole context
    // let (left_hyperkmer, right_hyperkmer) = extract_left_and_right_subsequences(
    //     hyperkmers,
    //     left_ext_hk_metadata,
    //     right_ext_hk_metadata,
    // );

    // let match_left = left_hk.equal_bitpacked(&left_hyperkmer);
    // let match_right = right_hk.equal_bitpacked(&right_hyperkmer);

    // match_left && match_right
}

/// Same as `search_exact_hyperkmer_match`, but checks the left first.
pub fn search_exact_hyperkmer_match_left_first(
    hyperkmers: &HyperkmerParts,
    left_hk: &Subsequence<NoBitPacked>,
    right_hk: &Subsequence<NoBitPacked>,
    left_ext_hk_metadata: &HKMetadata,
    right_ext_hk_metadata: &HKMetadata,
) -> bool {
    // extract left first
    let left_hyperkmer = extract_subsequence_from_hk_metadata(hyperkmers, left_ext_hk_metadata);
    // check left
    if !left_hk.equal_bitpacked(&left_hyperkmer) {
        return false;
    }

    // if matching on the left, only then check the right
    let right_hyperkmer = extract_subsequence_from_hk_metadata(hyperkmers, right_ext_hk_metadata);
    right_hk.equal_bitpacked(&right_hyperkmer)
}

/// Same as `search_exact_hyperkmer_match`, but checks the right first.
pub fn search_exact_hyperkmer_match_right_first(
    hyperkmers: &HyperkmerParts,
    left_hk: &Subsequence<NoBitPacked>,
    right_hk: &Subsequence<NoBitPacked>,
    left_ext_hk_metadata: &HKMetadata,
    right_ext_hk_metadata: &HKMetadata,
) -> bool {
    // extract right first
    let right_hyperkmer = extract_subsequence_from_hk_metadata(hyperkmers, right_ext_hk_metadata);
    // check right
    if !right_hk.equal_bitpacked(&right_hyperkmer) {
        return false;
    }

    // if matching on the right, only then check the left
    let left_hyperkmer = extract_subsequence_from_hk_metadata(hyperkmers, left_ext_hk_metadata);
    left_hk.equal_bitpacked(&left_hyperkmer)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyperkmer_contains_minimizer() {
        let minimizer: u64 = 56;
        let mut hk_count = HKCount::new();
        assert!(!hk_count.contains_minimizer(&minimizer));
        hk_count.insert_new_entry_in_hyperkmer_count(
            &minimizer,
            &HKMetadata::new(0, 0, 0, 5, true, false),
            &HKMetadata::new(0, 0, 0, 5, true, false),
            7,
        );
        hk_count.insert_new_entry_in_hyperkmer_count(
            &minimizer,
            &HKMetadata::new(0, 0, 0, 5, true, false),
            &HKMetadata::new(0, 0, 0, 5, true, false),
            7,
        );
        hk_count.insert_new_entry_in_hyperkmer_count(
            &(minimizer + 3),
            &HKMetadata::new(0, 0, 0, 5, true, false),
            &HKMetadata::new(0, 0, 0, 5, true, false),
            7,
        );
        hk_count.insert_new_entry_in_hyperkmer_count(
            &minimizer,
            &HKMetadata::new(0, 0, 0, 5, true, false),
            &HKMetadata::new(0, 0, 0, 5, true, false),
            7,
        );
        assert_eq!(
            hk_count.minimizer_set(),
            HashSet::from_iter(vec![minimizer, minimizer + 3])
        );

        assert!(!hk_count.get_data().is_empty());
    }
}
