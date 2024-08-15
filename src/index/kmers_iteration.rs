use crate::index::components::get_subsequence_from_metadata;
use crate::superkmer::SubsequenceMetadata;
use crate::Count;
use std::collections::HashMap;

use super::components::ExtendedHyperkmers;
use super::components::HKCount;
use super::components::HKMetadata;

pub fn extract_context_if_not_large(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
) -> Option<(String, usize)> {
    let (left_ext_hk_metadata, right_ext_hk_metadata, _count) = entry;
    if left_ext_hk_metadata.get_is_large() || right_ext_hk_metadata.get_is_large() {
        None
    } else {
        Some(extract_context(entry, m, hyperkmers, large_hyperkmers))
    }
}

pub fn extract_context_if_large(
    entry: &(HKMetadata, HKMetadata, Count),
    m: usize,
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
) -> Option<(String, usize)> {
    let (left_ext_hk_metadata, right_ext_hk_metadata, _count) = entry;
    if left_ext_hk_metadata.get_is_large() || right_ext_hk_metadata.get_is_large() {
        Some(extract_context(entry, m, hyperkmers, large_hyperkmers))
    } else {
        None
    }
}

pub fn extract_context_but_cut_the_minimizer_out(
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
    let left_context = left_ext_hk
        .subsequence(
            left_ext_hk_metadata.get_start(),
            left_ext_hk_metadata.get_end() - (m - 1),
        )
        .to_string();
    let right_context = right_ext_hk
        .subsequence(
            right_ext_hk_metadata.get_start() + (m - 1),
            right_ext_hk_metadata.get_end(),
        )
        .to_string();

    // let mut minimizer = left_ext_hk
    //     .subsequence(
    //         left_ext_hk_metadata.get_end() - (m - 1),
    //         left_ext_hk_metadata.get_end(),
    //     )
    //     .to_string();
    // minimizer.push_str(&right_ext_hk.subsequence(0, 1).to_string());

    // debug_assert_eq!(
    //     left_string[(left_string.len() - (m - 2))..left_string.len()],
    //     right_string[0..(m - 2)]
    // );

    let minimizer_pos = left_context.len();
    let whole_context = {
        let mut x = left_context;
        x.push_str(&right_context);
        x
    };

    (whole_context, minimizer_pos)
}

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

    let left_string = left_hyperkmer.to_string();
    let right_string = right_hyperkmer.to_string();

    debug_assert_eq!(
        left_string[(left_string.len() - (m - 2))..left_string.len()],
        right_string[0..(m - 2)]
    );
    let minimizer_pos = left_string.len() - (m - 1);
    let whole_context = {
        let mut x = left_string;
        x.push_str(&right_string[(m - 2)..right_string.len()]);
        x
    };

    // clearer with a name
    #[allow(clippy::let_and_return)]
    (whole_context, minimizer_pos)
}

fn extract_kmers_from_contexts_associated_to_a_minimizer(
    hk_count: &HKCount,
    minimizer: &u64,
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &[(usize, Vec<u8>)],
    k: usize,
    m: usize,
) -> Option<std::collections::hash_map::IntoIter<String, u16>> {
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
    // hk_count_iter: Box<dyn Iterator<Item = &'a (u64, (HKMetadata, HKMetadata, u16))> + 'a>,
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
        // let hk_count_iter = Box::new(hk_count.get_data().iter());

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
        // strategy: for each entry in the minimizer table for this minimizer, compute the original context
        // then extract kmers from this context
        // consider the union of the kmers from these contexts

        // let hk_count_elem = &self.hk_count_iter.next()?;
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

pub struct ContextsIterator<'a> {
    hk_count: &'a HKCount,
    minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
    // k: usize,
    m: usize,
    current_entry_iterator: Option<(
        u64,
        Box<dyn Iterator<Item = &'a (HKMetadata, HKMetadata, Count)> + 'a>,
    )>,
}

impl<'a> ContextsIterator<'a> {
    pub fn new(
        hk_count: &'a HKCount,
        minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
        hyperkmers: &'a ExtendedHyperkmers,
        large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
        // k: usize,
        m: usize,
    ) -> Self {
        // let hk_count_iter = Box::new(hk_count.get_data().iter());
        // let minimizer = minimizers_iter.peekable();

        Self {
            hk_count,
            minimizers_iter,
            hyperkmers,
            large_hyperkmers,
            // k,
            m,
            current_entry_iterator: None,
        }
    }
}

impl<'a> Iterator for ContextsIterator<'a> {
    type Item = (String, u64, usize, Count);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we have an iterator over entries, try to get the next item from it.
            if let Some((minimizer, entry_iter)) = &mut self.current_entry_iterator {
                if let Some(entry) = entry_iter.next() {
                    let context =
                        extract_context(entry, self.m, self.hyperkmers, self.large_hyperkmers);
                    return Some((context.0, *minimizer, context.1, entry.2));
                } else {
                    // If the current entry iterator is exhausted, clear it.
                    self.current_entry_iterator = None;
                }
            }

            // If there's no current entry iterator, try to create one from the next minimizer.
            if self.current_entry_iterator.is_none() {
                let minimizer = self.minimizers_iter.next()?;
                let hkcount_iter_for_this_minimizer = self.hk_count.get_data().get_iter(minimizer);
                // TODO eviter Box
                self.current_entry_iterator =
                    Some((*minimizer, Box::new(hkcount_iter_for_this_minimizer)));
            }
        }
    }
}

pub struct NormalContextsIterator<'a> {
    hk_count: &'a HKCount,
    minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
    // k: usize,
    m: usize,
    current_entry_iterator: Option<(
        u64,
        Box<dyn Iterator<Item = &'a (HKMetadata, HKMetadata, Count)> + 'a>,
    )>,
}

impl<'a> NormalContextsIterator<'a> {
    pub fn new(
        hk_count: &'a HKCount,
        minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
        hyperkmers: &'a ExtendedHyperkmers,
        large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
        // k: usize,
        m: usize,
    ) -> Self {
        // let hk_count_iter = Box::new(hk_count.get_data().iter());
        // let minimizer = minimizers_iter.peekable();

        Self {
            hk_count,
            minimizers_iter,
            hyperkmers,
            large_hyperkmers,
            // k,
            m,
            current_entry_iterator: None,
        }
    }
}

impl<'a> Iterator for NormalContextsIterator<'a> {
    type Item = (String, u64, usize, Count);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we have an iterator over entries, try to get the next item from it.
            if let Some((minimizer, entry_iter)) = &mut self.current_entry_iterator {
                if let Some(entry) = entry_iter.next() {
                    if let Some(context) = extract_context_if_not_large(
                        entry,
                        self.m,
                        self.hyperkmers,
                        self.large_hyperkmers,
                    ) {
                        return Some((context.0, *minimizer, context.1, entry.2));
                    } else {
                        continue;
                    }
                } else {
                    // If the current entry iterator is exhausted, clear it.
                    self.current_entry_iterator = None;
                }
            }

            // If there's no current entry iterator, try to create one from the next minimizer.
            if self.current_entry_iterator.is_none() {
                let minimizer = self.minimizers_iter.next()?;
                let hkcount_iter_for_this_minimizer = self.hk_count.get_data().get_iter(minimizer);
                // TODO eviter Box
                self.current_entry_iterator =
                    Some((*minimizer, Box::new(hkcount_iter_for_this_minimizer)));
            }
        }
    }
}

pub struct LargeContextsIterator<'a> {
    hk_count: &'a HKCount,
    minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
    m: usize,
    current_entry_iterator: Option<(
        u64,
        Box<dyn Iterator<Item = &'a (HKMetadata, HKMetadata, Count)> + 'a>,
    )>,
}

impl<'a> LargeContextsIterator<'a> {
    pub fn new(
        hk_count: &'a HKCount,
        minimizers_iter: std::collections::hash_set::IntoIter<&'a u64>,
        hyperkmers: &'a ExtendedHyperkmers,
        large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
        m: usize,
    ) -> Self {
        Self {
            hk_count,
            minimizers_iter,
            hyperkmers,
            large_hyperkmers,
            m,
            current_entry_iterator: None,
        }
    }
}

impl<'a> Iterator for LargeContextsIterator<'a> {
    type Item = (String, u64, usize, Count);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // If we have an iterator over entries, try to get the next item from it.
            if let Some((minimizer, entry_iter)) = &mut self.current_entry_iterator {
                if let Some(entry) = entry_iter.next() {
                    if let Some(context) = extract_context_if_large(
                        entry,
                        self.m,
                        self.hyperkmers,
                        self.large_hyperkmers,
                    ) {
                        return Some((context.0, *minimizer, context.1, entry.2));
                    } else {
                        continue;
                    }
                } else {
                    // If the current entry iterator is exhausted, clear it.
                    self.current_entry_iterator = None;
                }
            }

            // If there's no current entry iterator, try to create one from the next minimizer.
            if self.current_entry_iterator.is_none() {
                let minimizer = self.minimizers_iter.next()?;
                let hkcount_iter_for_this_minimizer = self.hk_count.get_data().get_iter(minimizer);
                // TODO eviter Box
                self.current_entry_iterator =
                    Some((*minimizer, Box::new(hkcount_iter_for_this_minimizer)));
            }
        }
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
