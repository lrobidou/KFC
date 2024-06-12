use crate::{
    get_rc_if_change_orientation, get_start_of_minimizer_in_superkmer,
    superkmers::compute_superkmers_linear, Count,
};

use super::{common_prefix_length, common_suffix_length, HKCount};

fn search_kmer(hk_count: &HKCount, hyperkmers: &[String], kmer: &str, k: usize, m: usize) -> Count {
    // TODO cost of doing k = kmer.len() ?
    let sks = compute_superkmers_linear(kmer, k, m);
    let superkmer = &sks[0]; // there can be only one
    let k = kmer.len();

    let match_size = k - 1;
    let mut total_count = 0;

    let start_of_minimizer_in_sk = get_start_of_minimizer_in_superkmer(superkmer);
    let left_sk = &superkmer.superkmer[0..(start_of_minimizer_in_sk + m - 1)];
    let right_sk = &superkmer.superkmer[(start_of_minimizer_in_sk + 1)..superkmer.superkmer.len()];

    for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count) in
        hk_count.get_iter(&superkmer.minimizer)
    {
        println!(
            "candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata = {:?}, {:?}",
            candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata
        );
        let candidate_left_ext_hk = get_rc_if_change_orientation(
            &hyperkmers[candidate_left_ext_hk_metadata.0],
            candidate_left_ext_hk_metadata.3,
        );
        let candidate_right_ext_hk = get_rc_if_change_orientation(
            &hyperkmers[candidate_right_ext_hk_metadata.0],
            candidate_right_ext_hk_metadata.3,
        );

        // extract candidate hyperkmers
        let candidate_left_hyperkmer = &candidate_left_ext_hk
            [candidate_left_ext_hk_metadata.1..candidate_left_ext_hk_metadata.2];
        let candidate_right_hyperkmer = &candidate_right_ext_hk
            [candidate_right_ext_hk_metadata.1..candidate_right_ext_hk_metadata.2];

        println!("candidate_left_ext_hk = {:?}", candidate_left_ext_hk);
        println!("candidate_right_ext_hk = {:?}", candidate_right_ext_hk);
        println!("candidate_left_hyperkmer = {:?}", candidate_left_hyperkmer);
        println!(
            "candidate_right_hyperkmer = {:?}",
            candidate_right_hyperkmer
        );
        let len_current_match_left = common_suffix_length(left_sk, candidate_left_hyperkmer);
        let len_current_match_right = common_prefix_length(right_sk, candidate_right_hyperkmer);
        let current_match_size = len_current_match_left + len_current_match_right;

        if current_match_size > match_size {
            total_count += count
        }
    }
    total_count
}

#[cfg(test)]
mod tests {
    // use crate::PosOfHyperkmerInExtHyperkmer;

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
        hk.insert(
            minimizer.into(),
            ((0, 0, 25, false), (1, 0, 31, false), count),
        );
        let search_result = search_kmer(&hk, &hyperkmers, kmer, kmer.len(), 10);
        assert!(search_result == count);
    }
}
