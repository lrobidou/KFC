use crate::superkmer::SubsequenceMetadata;
use std::collections::HashMap;

use super::super::components::HKCount;
use super::super::components::ParallelExtendedHyperkmers;
use super::extract_context;

pub fn extract_kmers_from_contexts_associated_to_a_minimizer(
    hk_count: &HKCount,
    minimizer: &u64,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
    k: &usize,
    m: &usize,
) -> std::collections::hash_map::IntoIter<Vec<u8>, u16> {
    // strategy: compute the original contexts associated with the minimizer
    // then extract kmers from these contexts
    // consider the union of the kmers from these contexts
    let mut kmers_counts: HashMap<Vec<u8>, u16> = HashMap::new();
    // for hk_count_chunk in hk_count.iter_chunks() {
    for hk_count_elem in hk_count.get_data().get_iter(minimizer) {
        let count = hk_count_elem.2;
        let (context, _minimizer_start_pos) =
            extract_context(hk_count_elem, *m, hyperkmers, large_hyperkmers);
        let context = SubsequenceMetadata::new(context.as_bytes(), 0, context.len(), true);
        // extract canonical kmers from the whole context and add them to the set
        for i in 0..context.len() - k + 1 {
            let kmer = context.subsequence(i, i + k).to_canonical().as_vec();

            kmers_counts
                .entry(kmer)
                .and_modify(|x| *x = x.saturating_add(count))
                .or_insert(count);
        }
    }
    kmers_counts.into_iter()
}
