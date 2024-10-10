mod hyperkmer_metadata;

use mashmap::MashMap;
use serde::{
    de::{MapAccess, Visitor},
    ser::SerializeMap,
    Deserialize, Deserializer, Serialize, Serializer,
};
use std::{
    collections::HashSet,
    sync::{Arc, RwLock},
};

use super::{add_new_hyperkmer, extended_hyperkmers::ParallelExtendedHyperkmers};
use super::{get_subsequence_from_metadata, ExtendedHyperkmers, LargeExtendedHyperkmer};
use crate::{
    check_equal_mashmap,
    subsequence::{BitPacked, NoBitPacked, Subsequence},
    Count, Minimizer,
};
use crate::{subsequence, Superkmer};

pub use hyperkmer_metadata::HKMetadata;

// Notation: in this module, variables starting with `c_` are `candidate`s variables,
// e.g. variables computed by iterating over the `HKCount` table.

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
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
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
            let hyperkmers = hyperkmers.get_bucket_from_id(bucket_id);
            let hyperkmers = hyperkmers.read().unwrap();
            let subseq = get_subsequence_from_metadata(
                &hyperkmers,
                large_hyperkmers,
                candidate_left_hk_metadata,
            );

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
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
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
            let hyperkmers = hyperkmers.get_bucket_from_id(bucket_id);
            let hyperkmers = hyperkmers.read().unwrap();

            let subseq = get_subsequence_from_metadata(
                &hyperkmers,
                large_hyperkmers,
                candidate_right_hk_metadata,
            );

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
        self.data
            .insert(*minimizer, (*left_metadata, *right_metadata, count));
    }

    /// Search if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`
    /// If so, increase the count of the occurence and return `true`
    /// Else, return false
    pub fn increase_count_if_exact_match(
        &mut self,
        minimizer: &Minimizer,
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
        left_hk: &Subsequence<NoBitPacked>,
        right_hk: &Subsequence<NoBitPacked>,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
            //DEBUG why not self mut ?
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

    /// Searches if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`.
    #[cfg(debug_assertions)]
    pub fn search_exact_match(
        &mut self,
        minimizer: &Minimizer,
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
        left_hk: &Subsequence<NoBitPacked>,
        right_hk: &Subsequence<NoBitPacked>,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
            self.data.get_iter(minimizer)
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
                return true;
            }
        }
        false
    }

    pub fn search_for_inclusion(
        &self,
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
        superkmer: &Superkmer,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let minimizer = superkmer.get_minimizer();
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count_hk) in
            self.data.get_iter(&minimizer)
        {
            // get sequences as they would appear if the current superkmer was canonical
            let is_large_left = candidate_left_ext_hk_metadata.get_is_large();
            let left_hyperkmers =
                hyperkmers.get_bucket_from_id(candidate_left_ext_hk_metadata.get_bucket_id());
            let left_hyperkmers = left_hyperkmers.read().unwrap();
            let subseq_left = get_subsequence_from_metadata(
                &left_hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            );
            let candidate_left_ext_hk: Subsequence<BitPacked> = subseq_left
                .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());

            let is_large_right = candidate_right_ext_hk_metadata.get_is_large();
            let right_hyperkmers =
                hyperkmers.get_bucket_from_id(candidate_right_ext_hk_metadata.get_bucket_id());
            let right_hyperkmers = right_hyperkmers.read().unwrap();
            let subseq_right = get_subsequence_from_metadata(
                &right_hyperkmers,
                large_hyperkmers,
                candidate_right_ext_hk_metadata,
            );
            let candidate_right_ext_hk: Subsequence<BitPacked> = subseq_right
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
                        candidate_left_ext_hk_metadata.get_bucket_id(),
                        candidate_left_ext_hk_metadata.get_index(),
                        start_left,
                        end_left,
                        is_large_left,
                        candidate_left_ext_hk_metadata.get_change_orientation(),
                    ),
                    HKMetadata::new(
                        candidate_right_ext_hk_metadata.get_bucket_id(),
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
    #[inline]
    pub fn contains_minimizer(&self, minimizer: &Minimizer) -> bool {
        self.data.contains_key(minimizer)
    }

    /// Search for the maximal inclusion of `left_sk` and `right_sk`
    /// in the extended hyperkmers associated with the minimizers
    /// Returns the metadata associated with this inclusion
    pub fn search_for_maximal_inclusion(
        &self,
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
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

        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let left_bucket_id = candidate_left_ext_hk_metadata.get_bucket_id();
            let left_hyperkmers = hyperkmers.get_bucket_from_id(left_bucket_id);
            let left_hyperkmers = left_hyperkmers.read().unwrap();
            let subseq_left = get_subsequence_from_metadata(
                &left_hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            );
            let candidate_left_ext_hk: Subsequence<BitPacked> = subseq_left
                .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());

            let right_bucket_id = candidate_right_ext_hk_metadata.get_bucket_id();
            let right_hyperkmers = hyperkmers.get_bucket_from_id(right_bucket_id);
            let right_hyperkmers = right_hyperkmers.read().unwrap();
            let subseq_right = get_subsequence_from_metadata(
                &right_hyperkmers,
                large_hyperkmers,
                candidate_right_ext_hk_metadata,
            );
            let candidate_right_ext_hk: Subsequence<BitPacked> = subseq_right
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

            debug_assert!(len_current_match_left >= m - 1);
            debug_assert!(len_current_match_right >= m - 1);

            if current_match_size - 2 * (m - 1) + m > match_size {
                match_size = current_match_size;

                match_metadata = Some((
                    // same suffix => same end, but different start
                    HKMetadata::new(
                        left_bucket_id,
                        candidate_left_ext_hk_metadata.get_index(),
                        candidate_left_ext_hk_metadata.get_end() - len_current_match_left,
                        candidate_left_ext_hk_metadata.get_end(),
                        candidate_left_ext_hk_metadata.get_is_large(),
                        candidate_left_ext_hk_metadata.get_change_orientation(),
                    ),
                    // same prefix => same start, but different end
                    HKMetadata::new(
                        right_bucket_id,
                        candidate_right_ext_hk_metadata.get_index(),
                        candidate_right_ext_hk_metadata.get_start(),
                        candidate_right_ext_hk_metadata.get_start() + len_current_match_right,
                        candidate_right_ext_hk_metadata.get_is_large(),
                        candidate_right_ext_hk_metadata.get_change_orientation(),
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
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
        minimizer: &Minimizer,
        left_context: &Subsequence<NoBitPacked>,
        right_context: &Subsequence<NoBitPacked>,
        k: usize,
        m: usize,
    ) -> Count {
        let mut total_count = 0;
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count) in
            self.data.get_iter(minimizer)
        {
            let left_hyperkmers_bucket_id = candidate_left_ext_hk_metadata.get_bucket_id();
            let left_hyperkmers = hyperkmers.get_bucket_from_id(left_hyperkmers_bucket_id);
            // TODO check for deadlocks
            let left_hyperkmers = left_hyperkmers.read().unwrap();
            // get sequences as they would appear if the current superkmer was canonical
            let candidate_left_hyperkmer = get_subsequence_from_metadata(
                &left_hyperkmers,
                large_hyperkmers,
                candidate_left_ext_hk_metadata,
            )
            .change_orientation_if(candidate_left_ext_hk_metadata.get_change_orientation());
            let candidate_left_hyperkmer = candidate_left_hyperkmer.subsequence(
                candidate_left_ext_hk_metadata.get_start(),
                candidate_left_ext_hk_metadata.get_end(),
            );

            let right_hyperkmers_bucket_id = candidate_right_ext_hk_metadata.get_bucket_id();
            let right_hyperkmers = hyperkmers.get_bucket_from_id(right_hyperkmers_bucket_id);
            let right_hyperkmers = right_hyperkmers.read().unwrap();
            let candidate_right_hyperkmer = get_subsequence_from_metadata(
                &right_hyperkmers,
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

    /// Finds a context that is an exact match on the right and an inclusion on the left. If this context exists, increments the count.
    pub fn increase_count_of_sk_if_found_else_insert_it(
        &mut self,
        k: usize,
        hyperkmers: &ParallelExtendedHyperkmers,
        large_hyperkmers: &[LargeExtendedHyperkmer],
        minimizer: &Minimizer,
        left_sk: &Subsequence<NoBitPacked>,
        right_sk: &Subsequence<NoBitPacked>,
    ) {
        let mut match_case = MatchCases::default();
        {
            // let mut found_match_right = false;
            for (c_left_hk_metadata, c_right_hk_metadata, count) in
                self.data.get_mut_iter(minimizer)
            {
                // extract ids
                let left_bucket_id = c_left_hk_metadata.get_bucket_id();
                let right_bucket_id = c_right_hk_metadata.get_bucket_id();

                // extract buckets
                let (left_hyperkmers, right_hyperkmers) =
                    hyperkmers.acquire_two_locks_read_mode(left_bucket_id, right_bucket_id);
                let right_hyperkmers = match right_hyperkmers.as_ref() {
                    Some(x) => x,
                    None => &left_hyperkmers,
                };

                // extract subsequences
                let (c_left_hk, c_right_hk) = extract_left_and_right_subsequences(
                    &left_hyperkmers,
                    right_hyperkmers,
                    large_hyperkmers,
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

                if match_case.is_min_cost() {
                    *count += 1;
                    break;
                }
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
    hyperkmers: &ParallelExtendedHyperkmers,
    k: usize,
) -> HKMetadata {
    // insert a new hyperkmer
    let nb_base_missing = (k - 1) - sequence.len();
    let metadata = if sequence.is_canonical() {
        assert!(sequence.is_canonical());
        let mut start = sequence.as_vec();
        start.extend_from_slice(&vec![b'A'; nb_base_missing]);
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        // TODO prove it
        assert!(full_sequence.is_canonical()); // Adding A at the end of a canonical read does not make it canonical

        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(&full_sequence);
        // create the associated metadata

        HKMetadata::new(id_bucket, id_hk, 0, sequence.len(), false, false)
    } else {
        let mut start = sequence.as_vec();
        start.extend_from_slice(&vec![b'T'; nb_base_missing]);
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        // TODO prove it
        assert!(!full_sequence.is_canonical()); // Adding T at the end of a non canonical read does not make it canonical
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
    hyperkmers: &ParallelExtendedHyperkmers,
    k: usize,
) -> HKMetadata {
    let nb_base_missing = (k - 1) - sequence.len();
    if sequence.is_canonical() {
        assert!(sequence.is_canonical());
        let mut start = vec![b'A'; nb_base_missing];
        start.extend_from_slice(&sequence.as_vec());
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        // TODO prove it
        assert!(full_sequence.is_canonical()); // Adding T at the beginning of a non canonical read does not make it canonical

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
        let full_sequence = Subsequence::new(&start, 0, k - 1, true);
        // TODO prove it
        assert!(!full_sequence.is_canonical()); // Adding T at the beginning of a non canonical read does not make it canonical
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

    fn is_min_cost(&self) -> bool {
        self.cost() == 0
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

pub fn extract_left_and_right_subsequences<'a>(
    left_hyperkmers: &'a ExtendedHyperkmers,
    right_hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a [LargeExtendedHyperkmer],
    left_hk_metadata: &HKMetadata,
    right_hk_metadata: &HKMetadata,
) -> (Subsequence<BitPacked<'a>>, Subsequence<BitPacked<'a>>) {
    // get sequences as they would appear if the current superkmer was canonical
    let left_hyperkmer =
        get_subsequence_from_metadata(left_hyperkmers, large_hyperkmers, left_hk_metadata)
            .change_orientation_if(left_hk_metadata.get_change_orientation());
    let left_hyperkmer =
        left_hyperkmer.subsequence(left_hk_metadata.get_start(), left_hk_metadata.get_end());

    let right_hyperkmer =
        get_subsequence_from_metadata(right_hyperkmers, large_hyperkmers, right_hk_metadata)
            .change_orientation_if(right_hk_metadata.get_change_orientation());
    let right_hyperkmer =
        right_hyperkmer.subsequence(right_hk_metadata.get_start(), right_hk_metadata.get_end());
    (left_hyperkmer, right_hyperkmer)
}

pub fn search_exact_hyperkmer_match(
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[LargeExtendedHyperkmer],
    left_hk: &Subsequence<NoBitPacked>,
    right_hk: &Subsequence<NoBitPacked>,
    left_ext_hk_metadata: &HKMetadata,
    right_ext_hk_metadata: &HKMetadata,
) -> bool {
    let left_bucket_id = left_ext_hk_metadata.get_bucket_id();
    let right_bucket_id = right_ext_hk_metadata.get_bucket_id();
    let (left_hyperkmers, right_hyperkmers) =
        hyperkmers.acquire_two_locks_read_mode(left_bucket_id, right_bucket_id);
    let right_hyperkmers = match right_hyperkmers.as_ref() {
        Some(x) => x,
        None => &left_hyperkmers,
    };

    // get sequences as they would appear if the current superkmer was canonical
    let left_hyperkmer =
        get_subsequence_from_metadata(&left_hyperkmers, large_hyperkmers, left_ext_hk_metadata);
    let left_hyperkmer =
        left_hyperkmer.change_orientation_if(left_ext_hk_metadata.get_change_orientation());
    let left_hyperkmer = left_hyperkmer.subsequence(
        left_ext_hk_metadata.get_start(),
        left_ext_hk_metadata.get_end(),
    );

    let right_hyperkmer =
        get_subsequence_from_metadata(right_hyperkmers, large_hyperkmers, right_ext_hk_metadata);
    let right_hyperkmer =
        right_hyperkmer.change_orientation_if(right_ext_hk_metadata.get_change_orientation());
    let right_hyperkmer = right_hyperkmer.subsequence(
        right_ext_hk_metadata.get_start(),
        right_ext_hk_metadata.get_end(),
    );

    let match_left = left_hk.equal_bitpacked(&left_hyperkmer);
    let match_right = right_hk.equal_bitpacked(&right_hyperkmer);

    match_left && match_right
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
