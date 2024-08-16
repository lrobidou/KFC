use crate::Count;

use super::super::components::ExtendedHyperkmers;
use super::super::components::HKCount;
use super::super::components::HKMetadata;
use super::extract_context;

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

pub struct ContextsIterator<'a> {
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

impl<'a> ContextsIterator<'a> {
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
