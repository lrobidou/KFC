use crate::index::components::extract_left_and_right_subsequences;
use crate::{index::components::get_subsequence_from_metadata, Count};

use super::components::{HKMetadata, LargeExtendedHyperkmer, ParallelExtendedHyperkmers};

use mashmap::IterGroupByKey;

use crate::subsequence::Subsequence;
use std::collections::HashMap;
use std::iter::Peekable;

pub fn extract_kmers_from_contexts_associated_to_a_minimizer(
    hk_counts_grouped_by_key: &mut Peekable<
        IterGroupByKey<'_, u64, (HKMetadata, HKMetadata, u16), ahash::RandomState>,
    >,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[LargeExtendedHyperkmer],
    k: &usize,
    m: &usize,
) -> Option<std::collections::HashMap<Vec<u8>, u16>> {
    let minimizer_for_this_loop = hk_counts_grouped_by_key.peek()?.0;
    // strategy: compute the original contexts associated with the minimizer
    // then extract kmers from these contexts
    // consider the union of the kmers from these contexts
    let mut kmers_counts: HashMap<Vec<u8>, u16> = HashMap::new();
    while let Some(next_entry) = hk_counts_grouped_by_key.peek() {
        // if the next minimizer is different than the current one, we return the current set
        if next_entry.0 != minimizer_for_this_loop {
            break;
        }
        let hk_count_elem = hk_counts_grouped_by_key.next().unwrap().1; // safe as peek was not None
        let count = hk_count_elem.2;
        let (context, _minimizer_start_pos) =
            extract_context(&hk_count_elem, *m, hyperkmers, large_hyperkmers);
        let context = Subsequence::new(context.as_bytes(), 0, context.len(), true);
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
    Some(kmers_counts)
}

pub fn extract_context(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[LargeExtendedHyperkmer],
) -> (String, usize) {
    let (left_ext_hk_metadata, right_ext_hk_metadata, _count) = entry;
    let (left_hyperkmers, right_hyperkmers) = hyperkmers.acquire_two_locks_read_mode(
        left_ext_hk_metadata.get_bucket_id(),
        right_ext_hk_metadata.get_bucket_id(),
    );
    let right_hyperkmers = match right_hyperkmers.as_ref() {
        Some(x) => x,
        None => &left_hyperkmers,
    };

    // extract relevant subsequences frome whole context
    let (left_hk, right_hk) = extract_left_and_right_subsequences(
        &left_hyperkmers,
        right_hyperkmers,
        large_hyperkmers,
        left_ext_hk_metadata,
        right_ext_hk_metadata,
    );

    let left_string = left_hk.as_vec();
    let right_string = right_hk.as_vec();

    // TODO there might be a way to prevent copy here
    debug_assert_eq!(
        left_string[(left_string.len() - (m - 2))..left_string.len()],
        right_string[0..(m - 2)]
    );
    let minimizer_pos = left_string.len() - (m - 1);
    let whole_context = {
        let mut x = left_string;
        let y = &right_string[(m - 2)..right_string.len()];
        x.extend_from_slice(y);
        x
    };

    (
        String::from_utf8(whole_context).expect("could not parse utf-8"),
        minimizer_pos,
    )
}
