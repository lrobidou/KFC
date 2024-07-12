use super::Superkmer;
use crate::{
    extended_hyperkmers::ExtendedHyperkmers,
    superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata},
    Count, Minimizer,
};
use mashmap::MashMap;

#[derive(Clone, Copy, Debug)]
pub struct HKMetadata {
    pub index: usize,
    pub start: usize,
    pub end: usize,
    // true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer it is associated with in the `HKCount` table
    pub change_orientation: bool,
}

pub struct HKCount {
    data: MashMap<Minimizer, (HKMetadata, HKMetadata, Count)>,
}

impl HKCount {
    pub fn new() -> Self {
        Self {
            data: MashMap::new(),
        }
    }

    /// Searches for a match between `extended_hyperkmer_left` and one of the left extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_left_id(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        minimizer: &Minimizer,
        extended_hyperkmer_left: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<usize> {
        let canonical_extended_hyperkmer_left = extended_hyperkmer_left.to_canonical();

        for (candidate_left_hk_metadata, _candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            if canonical_extended_hyperkmer_left.equal_bitpacked(
                &hyperkmers.get_hyperkmer_from_id(candidate_left_hk_metadata.index),
            ) {
                return Some(candidate_left_hk_metadata.index);
            }
        }
        None
    }

    /// Searches for a match between `extended_hyperkmer_right` and one of the right extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_right_id(
        &self,
        hyperkmers: &ExtendedHyperkmers,
        minimizer: &Minimizer,
        extended_hyperkmer_right: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<usize> {
        let canonical_extended_hyperkmer_right = extended_hyperkmer_right.to_canonical();
        for (_candidate_left_hk_metadata, candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            if canonical_extended_hyperkmer_right.equal_bitpacked(
                &hyperkmers.get_hyperkmer_from_id(candidate_right_hk_metadata.index),
            ) {
                return Some(candidate_right_hk_metadata.index);
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
        assert!(left_metadata.start < left_metadata.end);
        assert!(right_metadata.start < right_metadata.end);
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
        left_hk: &SubsequenceMetadata<NoBitPacked>,
        right_hk: &SubsequenceMetadata<NoBitPacked>,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
            self.data.get_mut_iter(minimizer)
        {
            let is_exact_match = search_exact_hyperkmer_match(
                hyperkmers,
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
        superkmer: &Superkmer,
        left_sk: &SubsequenceMetadata<NoBitPacked>,
        right_sk: &SubsequenceMetadata<NoBitPacked>,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let minimizer = superkmer.get_minimizer();
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count_hk) in
            self.data.get_mut_iter(&minimizer)
        {
            // get sequences as they would appear if the current superkmer was canonical
            let candidate_left_ext_hk: SubsequenceMetadata<BitPacked> = hyperkmers
                .get_hyperkmer_from_id(candidate_left_ext_hk_metadata.index)
                .change_orientation_if(candidate_left_ext_hk_metadata.change_orientation);

            let candidate_right_ext_hk: SubsequenceMetadata<BitPacked> = hyperkmers
                .get_hyperkmer_from_id(candidate_right_ext_hk_metadata.index)
                .change_orientation_if(candidate_right_ext_hk_metadata.change_orientation);

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
                    HKMetadata {
                        index: candidate_left_ext_hk_metadata.index,
                        start: start_left,
                        end: end_left,
                        change_orientation: candidate_left_ext_hk_metadata.change_orientation,
                    },
                    HKMetadata {
                        index: candidate_right_ext_hk_metadata.index,
                        start: start_right,
                        end: end_right,
                        change_orientation: candidate_right_ext_hk_metadata.change_orientation,
                    },
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
            // get sequences as they would appear if the current superkmer was canonical
            let candidate_left_ext_hk = hyperkmers
                .get_hyperkmer_from_id(candidate_left_ext_hk_metadata.index)
                .change_orientation_if(candidate_left_ext_hk_metadata.change_orientation);

            let candidate_right_ext_hk = hyperkmers
                .get_hyperkmer_from_id(candidate_right_ext_hk_metadata.index)
                .change_orientation_if(candidate_right_ext_hk_metadata.change_orientation);

            // extract candidate hyperkmers
            let candidate_left_hyperkmer = &candidate_left_ext_hk.subsequence(
                candidate_left_ext_hk_metadata.start,
                candidate_left_ext_hk_metadata.end,
            );
            let candidate_right_hyperkmer = &candidate_right_ext_hk.subsequence(
                candidate_right_ext_hk_metadata.start,
                candidate_right_ext_hk_metadata.end,
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
                    HKMetadata {
                        index: candidate_left_ext_hk_metadata.index,
                        start: candidate_left_ext_hk_metadata.end - len_current_match_left,
                        end: candidate_left_ext_hk_metadata.end,
                        change_orientation: candidate_left_ext_hk_metadata.change_orientation,
                    },
                    // same prefix => same start, but different end
                    HKMetadata {
                        index: candidate_right_ext_hk_metadata.index,
                        start: candidate_right_ext_hk_metadata.start,
                        end: candidate_right_ext_hk_metadata.start + len_current_match_right,
                        change_orientation: candidate_right_ext_hk_metadata.change_orientation,
                    },
                ));
            }
        }
        match_metadata
    }

    pub fn count_occurence_kmer(
        &self,
        hyperkmers: &ExtendedHyperkmers,
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
            let candidate_left_ext_hk = hyperkmers
                .get_hyperkmer_from_id(candidate_left_ext_hk_metadata.index)
                .change_orientation_if(candidate_left_ext_hk_metadata.change_orientation);

            let candidate_right_ext_hk = hyperkmers
                .get_hyperkmer_from_id(candidate_right_ext_hk_metadata.index)
                .change_orientation_if(candidate_right_ext_hk_metadata.change_orientation);

            // extract candidate hyperkmers
            let candidate_left_hyperkmer = &candidate_left_ext_hk.subsequence(
                candidate_left_ext_hk_metadata.start,
                candidate_left_ext_hk_metadata.end,
            );
            let candidate_right_hyperkmer = &candidate_right_ext_hk.subsequence(
                candidate_right_ext_hk_metadata.start,
                candidate_right_ext_hk_metadata.end,
            );
            let len_current_match_left =
                left_context.common_suffix_length_with_bitpacked(candidate_left_hyperkmer);
            let len_current_match_right =
                right_context.common_prefix_length_with_bitpacked(candidate_right_hyperkmer);
            let current_match_size = len_current_match_left + len_current_match_right;

            if current_match_size - 2 * (m - 1) + m >= k {
                total_count += count;
            }
        }
        total_count
    }
}

pub fn search_exact_hyperkmer_match(
    hyperkmers: &ExtendedHyperkmers,
    left_hk: &SubsequenceMetadata<NoBitPacked>,
    right_hk: &SubsequenceMetadata<NoBitPacked>,
    candidate_left_ext_hk_metadata: &HKMetadata,
    candidate_right_ext_hk_metadata: &HKMetadata,
) -> bool {
    // get sequences as they would appear if the current superkmer was canonical
    // TODO "style": this repetition of `candidate_left_hyperkmer` confuses me, how can I get rid of it ?
    let candidate_left_hyperkmer =
        hyperkmers.get_hyperkmer_from_id(candidate_left_ext_hk_metadata.index);
    let candidate_left_hyperkmer = candidate_left_hyperkmer
        .change_orientation_if(candidate_left_ext_hk_metadata.change_orientation);
    let candidate_left_hyperkmer = candidate_left_hyperkmer.subsequence(
        candidate_left_ext_hk_metadata.start,
        candidate_left_ext_hk_metadata.end,
    );

    let candidate_right_hyperkmer =
        hyperkmers.get_hyperkmer_from_id(candidate_right_ext_hk_metadata.index);
    let candidate_right_hyperkmer = candidate_right_hyperkmer
        .change_orientation_if(candidate_right_ext_hk_metadata.change_orientation);
    let candidate_right_hyperkmer = candidate_right_hyperkmer.subsequence(
        candidate_right_ext_hk_metadata.start,
        candidate_right_ext_hk_metadata.end,
    );

    let match_left = left_hk.equal_bitpacked(&candidate_left_hyperkmer);
    let match_right = right_hk.equal_bitpacked(&candidate_right_hyperkmer);

    match_left && match_right
}
