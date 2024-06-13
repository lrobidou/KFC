use crate::{
    common_prefix_length, common_suffix_length, get_rc_if_change_orientation,
    superkmers::{get_canonical_kmer, SuperKmerInfos},
    Count, Minimizer,
};
use mashmap::MashMap;

// (id of extended hk in the hyperkmers vector, start_pos, end_pos, orientation flag),
// orientation flag: true if the hyperkmer is in a DIFFERENT orientation
// that the canonical minimizer it is associated with in the `HKCount` table
type HKMetadata = (usize, usize, usize, bool);

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
        hyperkmers: &[String],
        minimizer: &str,
        extended_hyperkmer_left: &str,
    ) -> Option<usize> {
        let (_change_orientation, canonical_extended_hyperkmer_left) =
            get_canonical_kmer(extended_hyperkmer_left);
        for (candidate_left_hk_metadata, _candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            if canonical_extended_hyperkmer_left == hyperkmers[candidate_left_hk_metadata.0] {
                return Some(candidate_left_hk_metadata.0);
            }
        }
        None
    }

    /// Searches for a match between `extended_hyperkmer_right` and one of the right extended hyperkmers assiociated with `minimizer`
    /// minimizer is assumed to be in canonical form
    pub fn get_extended_hyperkmer_right_id(
        &self,
        hyperkmers: &[String],
        minimizer: &str,
        extended_hyperkmer_right: &str,
    ) -> Option<usize> {
        let (_change_orientation, canonical_extended_hyperkmer_right) =
            get_canonical_kmer(extended_hyperkmer_right);
        for (_candidate_left_hk_metadata, candidate_right_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            if canonical_extended_hyperkmer_right == hyperkmers[candidate_right_hk_metadata.0] {
                return Some(candidate_right_hk_metadata.0);
            }
        }
        None
    }

    /// Add a new entry in the hyperkmer counts.
    /// Minimizer should already be in canonical form.
    pub fn insert_new_entry_in_hyperkmer_count(
        &mut self,
        minimizer: &str,
        id_left_hk: (usize, usize, usize, bool),
        id_right_hk: (usize, usize, usize, bool),
        count: Count,
    ) {
        assert!(id_left_hk.1 < id_left_hk.2);
        assert!(id_right_hk.1 < id_right_hk.2);
        self.data
            .insert(String::from(minimizer), (id_left_hk, id_right_hk, count));
    }

    /// Search if `left_hk` and `right_hk` are associated with the minimizer of `superkmer`
    /// If so, increase the count of the occurence and return `true`
    /// Else, return false
    pub fn increase_count_if_exact_match(
        &self,
        superkmer: &SuperKmerInfos,
        hyperkmers: &[String],
        left_hk: &str,
        right_hk: &str,
    ) -> bool {
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
            self.data.get_mut_iter(&superkmer.minimizer)
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
        hyperkmers: &[String],
        superkmer: &SuperKmerInfos,
        left_sk: &str,
        right_sk: &str,
    ) -> Option<(HKMetadata, HKMetadata)> {
        // let mut new_left_and_right_metadata = None;
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count_hk) in
            self.data.get_mut_iter(&superkmer.minimizer)
        {
            // get sequences as they would appear if the current superkmer was canonical
            let candidate_left_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_left_ext_hk_metadata.0],
                candidate_left_ext_hk_metadata.3,
            );
            let candidate_right_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_right_ext_hk_metadata.0],
                candidate_right_ext_hk_metadata.3,
            );

            let match_start_left = candidate_left_ext_hk.starts_with(left_sk);
            let match_end_left = candidate_left_ext_hk.ends_with(left_sk);
            let match_left = match_start_left || match_end_left;

            let match_start_right = candidate_right_ext_hk.starts_with(right_sk);
            let match_end_right = candidate_right_ext_hk.ends_with(right_sk);
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
                    (
                        candidate_left_ext_hk_metadata.0,
                        start_left,
                        end_left,
                        candidate_left_ext_hk_metadata.3,
                    ),
                    (
                        candidate_right_ext_hk_metadata.0,
                        start_right,
                        end_right,
                        candidate_right_ext_hk_metadata.3,
                    ),
                ));
            }
        }
        None
    }

    /// Return true if the minimizer is present, false otherwise
    pub fn contains_minimizer(&self, minimizer: &str) -> bool {
        self.data.contains_key(minimizer)
    }

    /// Search for the maximal inclusion of `left_sk` and `right_sk`
    /// in the extended hyperkmers associated with the minimizers
    /// Returns the metadata associated with this inclusion
    pub fn search_for_maximal_inclusion(
        &self,
        hyperkmers: &[String],
        k: usize,
        minimizer: &str,
        left_sk: &str,
        right_sk: &str,
    ) -> Option<(HKMetadata, HKMetadata)> {
        let mut match_size = k - 1;
        let mut match_metadata = None;

        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
            self.data.get_iter(minimizer)
        {
            let candidate_left_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_left_ext_hk_metadata.0],
                candidate_left_ext_hk_metadata.3,
            );
            let candidate_right_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_right_ext_hk_metadata.0],
                candidate_right_ext_hk_metadata.3,
            );

            // extract candidate hyperkmers
            let candidate_left_hyperkmer = &candidate_left_ext_hk
                [candidate_left_ext_hk_metadata.1..candidate_left_ext_hk_metadata.2];
            let candidate_right_hyperkmer = &candidate_right_ext_hk
                [candidate_right_ext_hk_metadata.1..candidate_right_ext_hk_metadata.2];

            let len_current_match_left = common_suffix_length(left_sk, candidate_left_hyperkmer);
            let len_current_match_right = common_prefix_length(right_sk, candidate_right_hyperkmer);
            let current_match_size = len_current_match_left + len_current_match_right;

            assert!(len_current_match_left >= minimizer.len() - 1);
            assert!(len_current_match_right >= minimizer.len() - 1);

            if current_match_size - 2 * (minimizer.len() - 1) + minimizer.len() > match_size {
                match_size = current_match_size;

                match_metadata = Some((
                    // same suffix => same end, but different start
                    (
                        candidate_left_ext_hk_metadata.0,
                        candidate_left_ext_hk_metadata.2 - len_current_match_left,
                        candidate_left_ext_hk_metadata.2,
                        candidate_left_ext_hk_metadata.3,
                    ),
                    // same preffix => same start, but different end
                    (
                        candidate_right_ext_hk_metadata.0,
                        candidate_right_ext_hk_metadata.1,
                        candidate_right_ext_hk_metadata.1 + len_current_match_right,
                        candidate_right_ext_hk_metadata.3,
                    ),
                ));
            }
        }
        match_metadata
    }

    pub fn count_occurence_kmer(
        &self,
        hyperkmers: &[String],
        minimizer: &str,
        left_context: &str,
        right_context: &str,
        k: usize,
    ) -> Count {
        let mut total_count = 0;
        for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count) in
            self.data.get_iter(minimizer)
        {
            let candidate_left_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_left_ext_hk_metadata.0],
                candidate_left_ext_hk_metadata.3,
            );
            let candidate_right_ext_hk = get_rc_if_change_orientation(
                &hyperkmers[candidate_right_ext_hk_metadata.0],
                candidate_right_ext_hk_metadata.3,
            );

            // extract candidate hyperkmers
            let candidate_left_hyperkmer = &candidate_left_ext_hk
                [candidate_left_ext_hk_metadata.1..candidate_left_ext_hk_metadata.2];
            let candidate_right_hyperkmer = &candidate_right_ext_hk
                [candidate_right_ext_hk_metadata.1..candidate_right_ext_hk_metadata.2];
            let len_current_match_left =
                common_suffix_length(left_context, candidate_left_hyperkmer);
            let len_current_match_right =
                common_prefix_length(right_context, candidate_right_hyperkmer);
            let current_match_size = len_current_match_left + len_current_match_right;

            if current_match_size - 2 * (minimizer.len() - 1) + minimizer.len() >= k {
                total_count += count
            }
        }
        total_count
    }
}

fn search_exact_hyperkmer_match(
    hyperkmers: &[String],
    left_hk: &str,
    right_hk: &str,
    candidate_left_ext_hk_metadata: &(usize, usize, usize, bool),
    candidate_right_ext_hk_metadata: &(usize, usize, usize, bool),
) -> bool {
    // get sequences as they would appear if the current superkmer was canonical
    let candidate_left_ext_hk = get_rc_if_change_orientation(
        &hyperkmers[candidate_left_ext_hk_metadata.0],
        candidate_left_ext_hk_metadata.3,
    );
    let candidate_right_ext_hk = get_rc_if_change_orientation(
        &hyperkmers[candidate_right_ext_hk_metadata.0],
        candidate_right_ext_hk_metadata.3,
    );

    // extract candidate hyperkmers
    let candidate_left_hyperkmer =
        &candidate_left_ext_hk[candidate_left_ext_hk_metadata.1..candidate_left_ext_hk_metadata.2];
    let candidate_right_hyperkmer = &candidate_right_ext_hk
        [candidate_right_ext_hk_metadata.1..candidate_right_ext_hk_metadata.2];

    let match_left = candidate_left_hyperkmer == left_hk;
    let match_right = candidate_right_hyperkmer == right_hk;

    match_left && match_right
}
