use crate::{
    get_start_of_minimizer_in_superkmer, superkmers_computation::compute_superkmers_linear, Count,
};

use super::HKCount;

pub fn search_kmer(
    hk_count: &HKCount,
    hyperkmers: &[String],
    kmer: &str,
    k: usize,
    m: usize,
) -> Count {
    let sks = compute_superkmers_linear(kmer, k, m);
    let superkmer = &sks[0]; // there can be only one

    let start_of_minimizer_in_sk = get_start_of_minimizer_in_superkmer(superkmer);
    let left_sk = &superkmer.superkmer[0..(start_of_minimizer_in_sk + m - 1)];
    let right_sk = &superkmer.superkmer[(start_of_minimizer_in_sk + 1)..superkmer.superkmer.len()];

    hk_count.count_occurence_kmer(hyperkmers, &superkmer.minimizer, left_sk, right_sk, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search() {
        let kmer = "this_is_a_match_AAAAAAAAAA_this_is_another_match";
        let mut hk: HKCount = HKCount::new();
        let minimizer = "AAAAAAAAAA";
        let mut hyperkmers = Vec::new();
        let count = 34;

        let search_result = search_kmer(&hk, &hyperkmers, kmer, kmer.len(), 10);
        assert!(search_result == 0);

        hyperkmers.push(String::from("this_is_a_match_AAAAAAAAA"));
        hyperkmers.push(String::from("AAAAAAAAA_this_is_another_match"));
        hk.insert_new_entry_in_hyperkmer_count(
            minimizer,
            (0, 0, 25, false),
            (1, 0, 31, false),
            count,
        );
        let search_result = search_kmer(&hk, &hyperkmers, kmer, kmer.len(), 10);
        assert!(search_result == count);
    }
}
