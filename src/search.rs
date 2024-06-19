use crate::{
    compute_left_and_right::get_left_and_rigth_of_sk, superkmer::SubsequenceMetadata,
    superkmers_computation::compute_superkmers_linear, Count,
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
    let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

    hk_count.count_occurence_kmer(
        hyperkmers,
        &superkmer.get_minimizer(),
        &left_sk,
        &right_sk,
        k,
    )
}

#[cfg(test)]
mod tests {
    use crate::hyperkmers_counts::HKMetadata;

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
            &HKMetadata {
                index: 0,
                start: 0,
                end: 25,
                change_orientation: false,
            },
            &HKMetadata {
                index: 1,
                start: 0,
                end: 31,
                change_orientation: false,
            },
            count,
        );
        let search_result = search_kmer(&hk, &hyperkmers, kmer, kmer.len(), 10);
        assert!(search_result == count);
    }
}
