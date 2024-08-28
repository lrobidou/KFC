use mashmap::IterGroupByKey;

use crate::index::components::HKMetadata;
use crate::superkmer::SubsequenceMetadata;
use std::collections::HashMap;
use std::iter::Peekable;

use super::super::components::ParallelExtendedHyperkmers;
use super::extract_context;

pub fn extract_kmers_from_contexts_associated_to_a_minimizer(
    hk_counts_grouped_by_key: &mut Peekable<
        IterGroupByKey<'_, u64, (HKMetadata, HKMetadata, u16), ahash::RandomState>,
    >,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
    k: &usize,
    m: &usize,
) -> Option<std::collections::HashMap<Vec<u8>, u16>> {
    let minimizer_for_this_loop = hk_counts_grouped_by_key.peek()?.0;
    // strategy: compute the original contexts associated with the minimizer
    // then extract kmers from these contexts
    // consider the union of the kmers from these contexts
    let mut kmers_counts: HashMap<Vec<u8>, u16> = HashMap::new();
    // for hk_count_chunk in hk_count.iter_chunks() {
    while let Some(next_entry) = hk_counts_grouped_by_key.peek() {
        if next_entry.0 != minimizer_for_this_loop {
            return Some(kmers_counts);
        }
        let hk_count_elem = hk_counts_grouped_by_key.next().unwrap().1; // safe as peek was not None
        let count = hk_count_elem.2;
        let (context, _minimizer_start_pos) =
            extract_context(&hk_count_elem, *m, hyperkmers, large_hyperkmers);
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

    // end of the loop, nothing to see here
    None
}
