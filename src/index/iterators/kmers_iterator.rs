use crate::superkmer::SubsequenceMetadata;
use std::collections::HashMap;

use super::super::components::ExtendedHyperkmers;
use super::super::components::HKCount;
use super::extract_context;

fn extract_kmers_from_contexts_associated_to_a_minimizer(
    hk_count: &HKCount,
    minimizer: &u64,
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
    k: usize,
    m: usize,
) -> Option<std::collections::hash_map::IntoIter<String, u16>> {
    // strategy: compute the original contexts associated with the minimizer
    // then extract kmers from these contexts
    // consider the union of the kmers from these contexts
    let mut kmers_counts: HashMap<String, u16> = HashMap::new();
    for hk_count_elem in hk_count.get_data().get_iter(minimizer) {
        let count = hk_count_elem.2;
        let (context, _minimizer_start_pos) =
            extract_context(hk_count_elem, m, hyperkmers, large_hyperkmers);
        let context = SubsequenceMetadata::new(context.as_bytes(), 0, context.len(), true);
        // extract canonical kmers from the whole context and add them to the set
        for i in 0..context.len() - k + 1 {
            let kmer = context.subsequence(i, i + k).to_canonical().to_string();

            kmers_counts
                .entry(kmer)
                .and_modify(|x| *x = x.saturating_add(count))
                .or_insert(count);
        }
    }

    Some(kmers_counts.into_iter())
}

pub struct KmerIterator<'a> {
    hk_count: &'a HKCount,
    minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
    k: usize,
    m: usize,
}

impl<'a> KmerIterator<'a> {
    pub fn new(
        hk_count: &'a HKCount,
        minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
        hyperkmers: &'a ExtendedHyperkmers,
        large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
        k: usize,
        m: usize,
    ) -> Self {
        Self {
            hk_count,
            minimizers_iter,
            hyperkmers,
            large_hyperkmers,
            k,
            m,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = std::collections::hash_map::IntoIter<String, u16>;

    fn next(&mut self) -> Option<Self::Item> {
        let minimizer = self.minimizers_iter.next()?;
        extract_kmers_from_contexts_associated_to_a_minimizer(
            self.hk_count,
            minimizer,
            self.hyperkmers,
            self.large_hyperkmers,
            self.k,
            self.m,
        )
    }
}

// #[cfg(test)]
// mod tests {
// use itertools::Itertools;

// use crate::superkmers_computation;

// use super::*;

// #[test]
// fn test_extract_kmers_from_hk_count_iter_empty() {
//     let kmer = "TAAAGATCCCTCACAAAGAA";
//     let k = kmer.len();
//     let m = 10;
//     let superkmers =
//         superkmers_computation::compute_superkmers_linear_streaming(kmer.as_bytes(), k, m)
//             .expect("unable to extract superkmers")
//             .collect_vec();
//     assert_eq!(superkmers.len(), 1);
//     let superkmer = superkmers[0];
//     assert!(!superkmer.is_canonical_in_the_read());
//     let start_mini = superkmer.start_of_minimizer();

//     let left_context = &kmer[0..start_mini + m - 1];
//     let right_context = &kmer[start_mini + 1..kmer.len()];

//     let left_context =
//         SubsequenceMetadata::new(left_context.as_bytes(), 0, left_context.len(), false);
//     let right_context =
//         SubsequenceMetadata::new(right_context.as_bytes(), 0, right_context.len(), false);
// }
// }
