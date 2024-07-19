pub mod components;

use crate::Count;

use super::Minimizer;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use components::ExtendedHyperkmers;
use components::HKCount;
use components::SuperKmerCounts;

#[derive(Serialize, Deserialize, PartialEq)]
pub struct Index {
    super_kmer_counts: SuperKmerCounts,
    hk_count: HKCount,
    hyperkmers: ExtendedHyperkmers,
    discarded_minimizers: HashMap<Minimizer, u16>,
}

impl Index {
    pub fn new(
        super_kmer_counts: SuperKmerCounts,
        hk_count: HKCount,
        hyperkmers: ExtendedHyperkmers,
        discarded_minimizers: HashMap<u64, u16>,
    ) -> Self {
        Self {
            super_kmer_counts,
            hk_count,
            hyperkmers,
            discarded_minimizers,
        }
    }

    pub fn get_super_kmer_counts(&self) -> &SuperKmerCounts {
        &self.super_kmer_counts
    }
    pub fn get_hk_count(&self) -> &HKCount {
        &self.hk_count
    }
    pub fn get_hyperkmers(&self) -> &ExtendedHyperkmers {
        &self.hyperkmers
    }
    pub fn get_discarded_minimizers(&self) -> &HashMap<u64, u16> {
        &self.discarded_minimizers
    }

    pub fn search_kmer(&self, kmer: &[u8], k: usize, m: usize) -> Count {
        use crate::{
            compute_left_and_right::get_left_and_rigth_of_sk,
            superkmers_computation::compute_superkmers_linear_streaming,
        };
        let mut sks = compute_superkmers_linear_streaming(kmer, k, m).unwrap();
        let superkmer = &sks.next().unwrap(); // there can be only one
        debug_assert_eq!(sks.next(), None);
        let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

        self.hk_count.count_occurence_kmer(
            &self.hyperkmers,
            &superkmer.get_minimizer(),
            &left_sk,
            &right_sk,
            k,
            m,
        )
    }
}

#[cfg(test)]
mod tests {
    // use crate::{
    //     compute_left_and_right,
    //     hyperkmers_counts::HKMetadata,
    //     superkmer::SubsequenceMetadata, // two_bits::encode_minimizer,
    // };

    use crate::{
        compute_left_and_right, superkmer::SubsequenceMetadata,
        superkmers_computation::compute_superkmers_linear_streaming,
    };

    use super::*;
    use components::HKMetadata;
    use itertools::Itertools;

    #[test]
    fn test_search_empty_index() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 10;

        // nothing inserted => nothing is found
        let empty_index = Index::new(
            SuperKmerCounts::new(),
            HKCount::new(),
            ExtendedHyperkmers::new(k, 5),
            HashMap::new(),
        );
        let search_result = empty_index.search_kmer(kmer.as_bytes(), k, m);
        assert!(search_result == 0);
    }

    #[test]
    fn test_search() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();

        let superkmers = compute_superkmers_linear_streaming(kmer.as_bytes(), kmer.len(), 10)
            .unwrap()
            .collect_vec();
        assert!(superkmers.len() == 1);
        let superkmer = superkmers[0];

        // let's assume the minimizer is as follow
        assert_eq!(superkmer.start_of_minimizer(), 24);
        assert_eq!(superkmer.end_of_minimizer(), 34);

        let mut hk: HKCount = HKCount::new();
        let mut hyperkmers = ExtendedHyperkmers::new(kmer.len(), 7);
        let count = 34;

        // computing the left and right context to insert them in vector of hyperkmer
        let (left, right) = compute_left_and_right::get_left_and_rigth_of_sk(&superkmer);
        let size_left = left.len();
        let size_right = right.len();

        assert_eq!(size_left, superkmer.end_of_minimizer() - 1);
        assert_eq!(
            size_right,
            superkmer.superkmer.len() - (superkmer.start_of_minimizer() + 1)
        );

        // adding 'A' to complete the hyperkmers
        let mut left = left.to_string();
        while left.len() < k - 1 {
            left.push('A');
        }
        let mut right = right.to_string();
        while right.len() < k - 1 {
            right.push('A');
        }

        // computing the left and right extended hyperkmers
        let left = SubsequenceMetadata::new(left.as_bytes(), 0, left.len(), true);
        let right = SubsequenceMetadata::new(right.as_bytes(), 0, right.len(), true);

        // inserting the hyperkmers
        let index_left = hyperkmers.add_new_ext_hyperkmer(&left);
        let index_right = hyperkmers.add_new_ext_hyperkmer(&right);

        hk.insert_new_entry_in_hyperkmer_count(
            &superkmer.get_minimizer(),
            &HKMetadata {
                index: index_left,
                start: 0,
                end: size_left,
                change_orientation: left.is_canonical() != superkmer.is_canonical_in_the_read(),
            },
            &HKMetadata {
                index: index_right,
                start: 0,
                end: size_right,
                change_orientation: right.is_canonical() != superkmer.is_canonical_in_the_read(),
            },
            count,
        );

        let index = Index::new(SuperKmerCounts::new(), hk, hyperkmers, HashMap::new());
        let search_result = index.search_kmer(kmer.as_bytes(), kmer.len(), 10);
        assert_eq!(search_result, count);
    }
}
