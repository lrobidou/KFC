use crate::{index::components::get_subsequence_from_metadata, Count};

use super::components::{HKMetadata, ParallelExtendedHyperkmers};

use mashmap::IterGroupByKey;

use crate::subsequence::Subsequence;
use std::collections::HashMap;
use std::iter::Peekable;

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
    None
}

pub fn extract_context(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
) -> (String, usize) {
    let (left_ext_hk_metadata, right_ext_hk_metadata, _count) = entry;
    // get sequences as they would appear if the current superkmer was canonical
    let left_ext_hk =
        get_subsequence_from_metadata(hyperkmers, large_hyperkmers, left_ext_hk_metadata)
            .change_orientation_if(left_ext_hk_metadata.get_change_orientation());

    let right_ext_hk =
        get_subsequence_from_metadata(hyperkmers, large_hyperkmers, right_ext_hk_metadata)
            .change_orientation_if(right_ext_hk_metadata.get_change_orientation());

    // extract candidate hyperkmers
    let left_hyperkmer = &left_ext_hk.subsequence(
        left_ext_hk_metadata.get_start(),
        left_ext_hk_metadata.get_end(),
    );
    let right_hyperkmer = &right_ext_hk.subsequence(
        right_ext_hk_metadata.get_start(),
        right_ext_hk_metadata.get_end(),
    );

    let left_string = left_hyperkmer.as_vec();
    let right_string = right_hyperkmer.as_vec();

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

    // clearer with a name
    #[allow(clippy::let_and_return)]
    (
        String::from_utf8(whole_context).expect("could not parse utf-8"),
        minimizer_pos,
    )
}
