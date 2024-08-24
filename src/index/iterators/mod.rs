mod context_iterators;
mod kmers_iterator;

// pub use context_iterators::ContextsIterator;
pub use kmers_iterator::{extract_kmers_from_contexts_associated_to_a_minimizer, KmerIterator};

use crate::{index::components::get_subsequence_from_metadata, Count};

use super::components::{ExtendedHyperkmers, HKMetadata};

pub fn extract_context(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: &ExtendedHyperkmers,
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

fn guard_extract_context(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: std::sync::Arc<ExtendedHyperkmers>,
    large_hyperkmers: std::sync::Arc<[(usize, Vec<u8>)]>,
) -> (String, usize) {
    let (left_ext_hk_metadata, right_ext_hk_metadata, _count) = entry;
    // get sequences as they would appear if the current superkmer was canonical
    let left_ext_hk =
        get_subsequence_from_metadata(&hyperkmers, &large_hyperkmers, left_ext_hk_metadata)
            .change_orientation_if(left_ext_hk_metadata.get_change_orientation());

    let right_ext_hk =
        get_subsequence_from_metadata(&hyperkmers, &large_hyperkmers, right_ext_hk_metadata)
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
