use crate::index::components::get_subsequence_from_metadata;
use crate::superkmer::SubsequenceMetadata;
use std::collections::HashMap;
use std::collections::HashSet;

use super::components::ExtendedHyperkmers;
use super::components::HKCount;
use std::collections::hash_map;

fn extract_kmers_from_minimizer(
    hk_count: &HKCount,
    minimizer: &u64,
    hyperkmers: &ExtendedHyperkmers,
    large_hyperkmers: &Vec<(usize, Vec<u8>)>,
    k: usize,
    m: usize,
) -> Option<hash_map::IntoIter<String, u16>> {
    // strategy: for each entry in the minimizer table for this minimizer, compute the original context
    // then extract kmers from this context
    // consider the union of the kmers from these contexts
    let mut kmers_counts: HashMap<String, u16> = HashMap::new();
    for hk_count_elem in hk_count.get_data().get_iter(minimizer) {
        let (left_ext_hk_metadata, right_ext_hk_metadata, count) = hk_count_elem;
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

        let whole_context = {
            let mut x = left_string;
            x.push_str(&right_string[(m - 2)..right_string.len()]);
            x
        };

        let whole_context =
            SubsequenceMetadata::new(whole_context.as_bytes(), 0, whole_context.len(), true);

        // extract canonical kmers from the whole context and add them to the set
        for i in 0..whole_context.len() - k + 1 {
            let kmer = whole_context.subsequence(i, i + k);
            let kmer = kmer.to_canonical();
            let kmer = kmer.to_string();

            kmers_counts
                .entry(kmer)
                .and_modify(|x| *x = x.saturating_add(*count))
                .or_insert(*count);
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
        // minimizer_iter: std::collections::hash_set::Iter<&u64>,
        hyperkmers: &'a ExtendedHyperkmers,
        large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
        k: usize,
        m: usize,
    ) -> Self {
        println!("new kmers iterator");

        let mut set = HashSet::new();
        for (k, _v) in hk_count.get_data().iter() {
            set.insert(k);
        }
        let minimizers_iter = set.into_iter();

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
        // let hk_count_elem = &self.hk_count_iter.next()?;
        let minimizer = self.minimizers_iter.next()?;
        extract_kmers_from_minimizer(
            self.hk_count,
            minimizer,
            self.hyperkmers,
            self.large_hyperkmers,
            self.k,
            self.m,
        )
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::superkmers_computation;

    use super::*;

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
}
