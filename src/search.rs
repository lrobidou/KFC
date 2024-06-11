use crate::{
    get_rc_if_change_orientation,
    superkmers::{compute_superkmers_linear, reverse_complement},
    Count,
};

use super::{
    common_prefix_length, common_suffix_length, get_rc_if_change_orientation_internal, HKCount,
};

fn search_kmer(hk: &HKCount, hyperkmers: &[String], kmer: &str, k: usize, m: usize) -> Count {
    // TODO cost of doing k = kmer.len() ?
    let sks = compute_superkmers_linear(kmer, k, m);
    let sk = &sks[0]; // there can be only one

    let start_of_minimizer_in_sk = sk.start_of_minimizer_as_read - sk.start_of_superkmer_as_read;
    let left_hk_of_sk = &kmer[0..(start_of_minimizer_in_sk + m - 1)];
    let right_hk_of_sk = &kmer[(start_of_minimizer_in_sk + 1)..kmer.len()];
    // println!("{:?}", sk);
    // println!("searching for {left_hk_of_sk}, {right_hk_of_sk}");
    for (
        (left_id, left_size, left_orientation),
        (right_id, right_size, right_orientation),
        count,
    ) in hk.get_iter(&sk.minimizer)
    {
        let candidate_left_hk =
            get_rc_if_change_orientation(&hyperkmers[*left_id], *left_orientation);
        let candidate_right_hk =
            get_rc_if_change_orientation(&hyperkmers[*right_id], *right_orientation);

        let is_match_left =
            common_suffix_length(left_hk_of_sk, &candidate_left_hk) == left_hk_of_sk.len();
        let is_match_right =
            common_prefix_length(right_hk_of_sk, &candidate_right_hk) == right_hk_of_sk.len();

        // println!("{:?}", (&candidate_left_hk, &candidate_right_hk));
        // println!("{:?}", (&is_match_left, &is_match_right));

        if is_match_left && is_match_right {
            return *count;
        }
    }
    0
}

// // BUG when I search on the left, I should also check the rigth if the minimizer is not canonical

// /// search a sequence at the left of the minimizer in the index.
// /// Minimizer and sequence must be in canonical form.
// /// If minimizer is not found, then returns false.
// fn search_left_inner<ReversComplementFunction>(
//     hk: &HKCount,
//     hyperkmers: &[String],
//     m: usize,
//     minimizer: &str,
//     seq: &str,
//     minimizer_start_pos: usize,
//     right_minimizer_start_pos: usize,
//     revcompfunc: &ReversComplementFunction,
// ) -> bool
// where
//     ReversComplementFunction: Fn(&str) -> String,
// {
//     let left_context = &seq[0..minimizer_start_pos + m - 1]; // TODO check
//     let right_context = &seq[minimizer_start_pos + 1..right_minimizer_start_pos - 1];
//     // we start by ensuring we only select tuple (left, right) where we match the right
//     // TODO I did thuis snippet way too often.
//     // I should provide a way to iterate on hyperkmers with the correct orientation
//     for ((id_left, change_orientation_left), (id_right, change_orientation_right), _count) in
//         hk.get_iter(minimizer)
//     {
//         let candidate_left_hk = get_rc_if_change_orientation_internal(
//             revcompfunc,
//             &hyperkmers[*id_left],
//             *change_orientation_left,
//         );
//         let candidate_right_hk = get_rc_if_change_orientation_internal(
//             revcompfunc,
//             &hyperkmers[*id_right],
//             *change_orientation_right,
//         );

//         println!("{:?}", (&left_context, &right_context));
//         println!("{:?}", (&candidate_left_hk, &candidate_right_hk));

//         // check if the sequence at the right of the minimizer is matching the right hyperkemr candidate
//         let len_match_right = common_prefix_length(right_context, &candidate_right_hk);
//         if len_match_right != std::cmp::min(right_context.len(), candidate_right_hk.len()) {
//             // mismatch, see you next iteration
//             println!("continuing {len_match_right} != min({right_context}, {candidate_right_hk})");
//             continue;
//         }

//         let len_match_left = common_suffix_length(left_context, &candidate_left_hk);
//         if len_match_left == left_context.len() {
//             // our entire context is found for a minimizer that matches on the right
//             return true;
//         } else if len_match_left == candidate_left_hk.len() {
//             // we stumble accros a minimizer
//             let new_minimizer_start_pos = minimizer_start_pos - len_match_left - 1;
//             let new_minimizer = &seq[new_minimizer_start_pos..new_minimizer_start_pos + m];

//             if search_left_inner(
//                 hk,
//                 hyperkmers,
//                 m,
//                 new_minimizer,
//                 seq,
//                 new_minimizer_start_pos,
//                 minimizer_start_pos,
//                 revcompfunc,
//             ) {
//                 return true;
//             } else {
//                 // mismatch
//                 // maybe we'll have more chance at the next iteration ?
//                 continue;
//             }
//         } else {
//             // mismatch
//             // maybe we'll have more chance at the next iteration ?
//             continue;
//         }
//     }
//     false
// }

// /// search a sequence at the right of the minimizer in the index.
// /// Minimizer and sequence must be in canonical form.
// /// If minimizer is not found, then returns false.
// fn search_right<ReversComplementFunction>(
//     hk: &HKCount,
//     hyperkmers: &[String],
//     m: usize,
//     minimizer: &str,
//     seq: &str,
//     minimizer_start_pos: usize,
//     left_minimizer_start_pos: usize,
//     revcompfunc: &ReversComplementFunction,
// ) -> bool
// where
//     ReversComplementFunction: Fn(&str) -> String,
// {
//     let left_context = &seq[left_minimizer_start_pos + 1..minimizer_start_pos + m - 1];
//     let right_context = &seq[minimizer_start_pos + 1..seq.len()];
//     // we start by ensuring we only select tuple (left, right) where we match the left
//     // TODO I did this snippet way too often.
//     // I should provide a way to iterate on hyperkmers with the correct orientation
//     for ((id_left, change_orientation_left), (id_right, change_orientation_right), _count) in
//         hk.get_iter(minimizer)
//     {
//         // OPTIMIZE verbeux
//         let candidate_left_hk = if *change_orientation_left {
//             revcompfunc(&hyperkmers[*id_left])
//         } else {
//             String::from(&hyperkmers[*id_left])
//         };
//         let candidate_right_hk = if *change_orientation_right {
//             revcompfunc(&hyperkmers[*id_right])
//         } else {
//             String::from(&hyperkmers[*id_right])
//         };

//         // check if the sequence at the left of the minimizer is matching the left hyperkemr candidate
//         let len_match_left = common_suffix_length(left_context, &candidate_left_hk);
//         if len_match_left != std::cmp::min(left_context.len(), candidate_left_hk.len()) {
//             // mismatch, see you next iteration
//             continue;
//         }

//         let len_match_right = common_prefix_length(right_context, &candidate_right_hk);
//         if len_match_right == right_context.len() {
//             // our entire context is found for a minimizer that matches on the left
//             return true;
//         } else if len_match_right == candidate_right_hk.len() {
//             // we stumble accros a minimizer
//             let new_minimizer_start_pos = minimizer_start_pos + (len_match_right + 1) + 1;
//             let new_minimizer = &seq[new_minimizer_start_pos..new_minimizer_start_pos + m];

//             if search_right(
//                 hk,
//                 hyperkmers,
//                 m,
//                 new_minimizer,
//                 seq,
//                 new_minimizer_start_pos,
//                 minimizer_start_pos,
//                 revcompfunc,
//             ) {
//                 return true;
//             } else {
//                 // mismatch
//                 // maybe we'll have more chance at the next iteration ?
//                 continue;
//             }
//         } else {
//             // mismatch
//             // maybe we'll have more chance at the next iteration ?
//             continue;
//         }
//     }
//     false
// }

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
        hk.insert(minimizer.into(), ((0, 14, false), (1, 30, false), count));
        let search_result = search_kmer(&hk, &hyperkmers, kmer, kmer.len(), 10);
        assert!(search_result == count);
    }
}
