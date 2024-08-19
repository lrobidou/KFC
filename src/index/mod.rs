pub mod components;
mod computation;
mod iterators;

use super::Minimizer;
use crate::Count;
use computation::first_stage;
use computation::second_stage;
use iterators::ContextsIterator;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;
use std::time::Instant;

use iterators::KmerIterator;

use components::ExtendedHyperkmers;
use components::HKCount;
use components::SuperKmerCounts;

use crate::{
    compute_left_and_right::get_left_and_rigth_of_sk,
    superkmers_computation::compute_superkmers_linear_streaming,
};
// -- type state pattern
pub trait FullIndexTrait {}

#[derive(Serialize, Deserialize, PartialEq)]
pub struct CompleteIndex {
    super_kmer_counts: SuperKmerCounts,
    discarded_minimizers: HashMap<Minimizer, u16>,
}

#[derive(Serialize, Deserialize, PartialEq)]
pub struct StrippedIndex {}

impl FullIndexTrait for CompleteIndex {}
impl FullIndexTrait for StrippedIndex {}

/// Index of KFC
/// The index holds the hyperkmers and their counts.
/// During contruction, the index also holds the counts of superkmers and information about minimizers.
#[derive(Serialize, Deserialize, PartialEq)]
pub struct Index<FI: FullIndexTrait> {
    hk_count: HKCount,
    hyperkmers: ExtendedHyperkmers,
    /// vector of larger extended hyperkmers // TODO document
    large_hyperkmers: Vec<(usize, Vec<u8>)>, // TODO use slice
    k: usize,
    m: usize,
    superkmers_infos: FI,
}

impl Index<CompleteIndex> {
    /// Constructs a new `Index` indexing a set of sequences
    #[allow(clippy::self_named_constructors)] // Self named constructor ? I want it that way 🎵
    pub fn index(k: usize, m: usize, threshold: Count, sequences: &Vec<&str>) -> Self {
        let start_fisrt_step = Instant::now();
        let (mut super_kmer_counts, mut hk_count, hyperkmers, large_hyperkmers) =
            first_stage(sequences, k, m, threshold);
        println!(
            "time first stage: {} milliseconds",
            start_fisrt_step.elapsed().as_millis()
        );
        let start_second_stage = Instant::now();
        let discarded_minimizers = second_stage(
            &mut super_kmer_counts,
            &mut hk_count,
            &hyperkmers,
            &large_hyperkmers,
            sequences,
            k,
            m,
            threshold,
        );
        println!(
            "time second stage: {} milliseconds",
            start_second_stage.elapsed().as_millis()
        );
        Self {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            superkmers_infos: CompleteIndex {
                super_kmer_counts,
                discarded_minimizers,
            },
            k,
            m,
        }
    }

    pub fn remove_superkmer_infos(self) -> Index<StrippedIndex> {
        Index {
            hk_count: self.hk_count,
            hyperkmers: self.hyperkmers,
            large_hyperkmers: self.large_hyperkmers,
            k: self.k,
            m: self.m,
            superkmers_infos: StrippedIndex {},
        }
    }
}

impl Index<StrippedIndex> {
    /// Constructs a new `Index` indexing a set of sequences
    #[allow(clippy::self_named_constructors)] // Self named constructor ? I want it that way 🎵
    pub fn index(k: usize, m: usize, threshold: Count, sequences: &Vec<&str>) -> Self {
        let start_fisrt_step = Instant::now();
        let (mut super_kmer_counts, mut hk_count, hyperkmers, large_hyperkmers) =
            first_stage(sequences, k, m, threshold);
        println!(
            "time first stage: {} milliseconds",
            start_fisrt_step.elapsed().as_millis()
        );
        let start_second_stage = Instant::now();
        let _discarded_minimizers = second_stage(
            &mut super_kmer_counts,
            &mut hk_count,
            &hyperkmers,
            &large_hyperkmers,
            sequences,
            k,
            m,
            threshold,
        );
        println!(
            "time second stage: {} milliseconds",
            start_second_stage.elapsed().as_millis()
        );
        Self {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            superkmers_infos: StrippedIndex {},
            k,
            m,
        }
    }
}

impl<FI> Index<FI>
where
    FI: FullIndexTrait + Serialize,
{
    #[cfg(test)]
    pub fn new(
        hk_count: HKCount,
        hyperkmers: ExtendedHyperkmers,
        large_hyperkmers: Vec<(usize, Vec<u8>)>,
        superkmers_infos: FI,
        k: usize,
        m: usize,
    ) -> Self {
        Self {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            superkmers_infos,
            k,
            m,
        }
    }

    pub fn get_hyperkmers(&self) -> &ExtendedHyperkmers {
        &self.hyperkmers
    }

    pub fn get_large_hyperkmers(&self) -> &Vec<(usize, Vec<u8>)> {
        &self.large_hyperkmers
    }

    pub fn search_kmer(&self, kmer: &[u8], k: usize, m: usize) -> Count {
        // there can be only one superkmer
        let mut sks = compute_superkmers_linear_streaming(kmer, k, m).unwrap();
        let superkmer = &sks.next().unwrap();
        debug_assert_eq!(sks.next(), None);

        let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

        self.hk_count.count_occurence_kmer(
            &self.hyperkmers,
            &self.large_hyperkmers,
            &superkmer.get_minimizer(),
            &left_sk,
            &right_sk,
            k,
            m,
        )
    }

    pub fn iter_kmers(&self) -> KmerIterator {
        let minimizers = self.iter_minimizers();
        KmerIterator::new(
            &self.hk_count,
            minimizers,
            &self.hyperkmers,
            &self.large_hyperkmers,
            self.k,
            self.m,
        )
    }

    pub fn iter_minimizers(&self) -> std::collections::hash_set::IntoIter<&u64> {
        let mut set = HashSet::new();
        for (minimizer, _v) in self.hk_count.get_data().iter() {
            set.insert(minimizer);
        }
        set.into_iter()
    }

    /// Once a change in minimizer occurs, it never comes back again
    pub fn context_iterator(&self) -> ContextsIterator {
        let mini_iter = self.iter_minimizers();
        ContextsIterator::new(
            &self.hk_count,
            mini_iter,
            &self.hyperkmers,
            &self.large_hyperkmers,
            self.m,
        )
    }

    // pub fn iter_minimizers_large_context(&self) -> std::collections::hash_set::IntoIter<&u64> {
    //     let mut set = HashSet::new();
    //     for (minimizer, (l, r, _c)) in self.hk_count.get_data().iter() {
    //         if l.get_is_large() || r.get_is_large() {
    //             set.insert(minimizer);
    //         }
    //     }
    //     set.into_iter()
    // }

    // TODO
    // pub fn iter_minimizers_normal_contexts(&self) -> std::collections::hash_set::IntoIter<&u64> {
    //     let mut set = HashSet::new();
    //     for (k, _v) in self.hk_count.get_data().iter() {
    //         set.insert(k);
    //     }
    //     set.into_iter()
    // }
    pub fn get_k(&self) -> usize {
        self.k
    }

    pub fn get_m(&self) -> usize {
        self.m
    }
}

#[cfg(test)]
mod tests {

    use crate::serde::bin;

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
        let empty_index: Index<CompleteIndex> = Index::new(
            HKCount::new(),
            ExtendedHyperkmers::new(k, 5),
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: SuperKmerCounts::new(),
                discarded_minimizers: HashMap::new(),
            },
            k,
            m,
        );
        let search_result = empty_index.search_kmer(kmer.as_bytes(), k, m);
        assert!(search_result == 0);
    }

    #[test]
    fn test_search() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 10;

        let superkmers = compute_superkmers_linear_streaming(kmer.as_bytes(), kmer.len(), m)
            .unwrap()
            .collect_vec();
        assert!(superkmers.len() == 1);
        let superkmer = superkmers[0];

        // let's assume the minimizer is as follow
        assert_eq!(superkmer.start_of_minimizer(), 24);
        assert_eq!(superkmer.end_of_minimizer(), 34);

        let mut hk: HKCount = HKCount::new();
        let mut hyperkmers = ExtendedHyperkmers::new(kmer.len(), 7);
        let large_hyperkmers = Vec::new();
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

        let left_hk_metadata = HKMetadata::new(
            index_left,
            0,
            size_left,
            false,
            left.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let right_hk_metadata = HKMetadata::new(
            index_right,
            0,
            size_right,
            false,
            right.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        hk.insert_new_entry_in_hyperkmer_count(
            &superkmer.get_minimizer(),
            &left_hk_metadata,
            &right_hk_metadata,
            count,
        );

        let index: Index<CompleteIndex> = Index::new(
            hk,
            hyperkmers,
            large_hyperkmers,
            CompleteIndex {
                super_kmer_counts: SuperKmerCounts::new(),
                discarded_minimizers: HashMap::new(),
            },
            k,
            m,
        );
        let search_result = index.search_kmer(kmer.as_bytes(), kmer.len(), 10);
        assert_eq!(search_result, count);
    }

    #[test]
    fn test_empty_serde() {
        let filename = "test_empty_serde";
        let k = 10;
        let m = 3;

        let empty_index: Index<CompleteIndex> = Index::new(
            HKCount::new(),
            ExtendedHyperkmers::new(k, 5),
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: SuperKmerCounts::new(),
                discarded_minimizers: HashMap::new(),
            },
            k,
            m,
        );
        let empty_index = empty_index.remove_superkmer_infos();
        bin::dump(&empty_index, filename).unwrap();
        assert!(empty_index == bin::load(filename).unwrap());
    }

    #[test]
    fn test_serde() {
        let filename = "test_serde";
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 3;

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

        let left_hk_metadata = HKMetadata::new(
            index_left,
            0,
            size_left,
            false,
            left.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let right_hk_metadata = HKMetadata::new(
            index_right,
            0,
            size_right,
            false,
            right.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        hk.insert_new_entry_in_hyperkmer_count(
            &superkmer.get_minimizer(),
            &left_hk_metadata,
            &right_hk_metadata,
            count,
        );

        let index: Index<CompleteIndex> = Index::new(
            hk,
            hyperkmers,
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: SuperKmerCounts::new(),
                discarded_minimizers: HashMap::new(),
            },
            k,
            m,
        );

        let index = index.remove_superkmer_infos();
        bin::dump(&index, filename).unwrap();
        assert!(index == bin::load(filename).unwrap());
    }

    #[test]
    fn test_iter_kmers_empty() {
        let index: Index<CompleteIndex> = Index::<CompleteIndex>::index(21, 11, 1, &vec![]);
        let kmers_iter = index.iter_kmers().flatten();
        let kmers = kmers_iter.collect_vec();
        assert_eq!(kmers, vec![]);
    }

    // #[test]
    // fn test_iter_kmers_one_element() {
    //     let v: HashSet<(String, Count)> = HashSet::from_iter(
    //         vec![
    //             ("TAGCTTCTCGCTATTAGCTTC", 1),
    //             ("AGCTTCTCGCTATTAGCTTCA", 1),
    //             ("GCTTCTCGCTATTAGCTTCAA", 1),
    //             ("CTTCTCGCTATTAGCTTCAAT", 1),
    //             ("TTCTCGCTATTAGCTTCAATG", 1),
    //             ("TCTCGCTATTAGCTTCAATGA", 1),
    //             ("CTCGCTATTAGCTTCAATGAT", 1),
    //             ("TCGCTATTAGCTTCAATGATA", 1),
    //             ("CGCTATTAGCTTCAATGATAC", 1),
    //             ("GCTATTAGCTTCAATGATACG", 1),
    //         ]
    //         .into_iter()
    //         .map(|x| (x.0.into(), x.1)),
    //     );

    //     let index = Index::index(21, 11, 1, &vec!["TAGCTTCTCGCTATTAGCTTCAATGATACG"]);
    //     let kmers_results: HashSet<(String, Count)> = HashSet::from_iter(index.iter_kmers());
    //     assert_eq!(kmers_results, v);
    // }
}
