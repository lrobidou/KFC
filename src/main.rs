use compute_left_and_right::{get_left_and_rigth_extended_hk, get_left_and_rigth_hk};
use fastxgz::fasta_reads;
use log::warn;
use mashmap::MashMap;
use std::collections::HashMap;
use std::time::Instant;

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;
// type SKCount = MashMap<Minimizer, (HashSuperKmer, Count)>;

#[derive(Debug, Clone, Copy)]
enum PosOfHyperkmerInExtHyperkmer {
    Start,
    End,
}

// canonical minimmizer -> (
//     (id of left extended hk, hk len, StartingPos, orientation flag),
//     (id of right extended hk, hk len, StartingPos, orientation flag),
//     count
// )
// orientation flag: true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer
type HKCount = MashMap<
    Minimizer,
    (
        (usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
        (usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
        Count,
    ),
>;
// same as HKCount
// but position is now indicated by (id, start_pos, end_pos, orientation flag)
// (i.e. no more len and StartingPos)
type TruncatedHKCount = MashMap<
    Minimizer,
    (
        (usize, usize, usize, bool),
        (usize, usize, usize, bool),
        Count,
    ),
>;

mod brrr_minimizers;
mod compute_left_and_right;
mod search;
mod superkmers;

use superkmers::{
    compute_superkmers_linear, get_canonical_kmer, reverse_complement, SuperKmerInfos,
};

mod superkmers_count;
use superkmers_count::SuperKmerCounts;

/// Searches for `extended_hyperkmer_left` in `hk_count[minimizer]`
/// minimizer is assumed to be in canonical form
fn get_extended_hyperkmer_left_id(
    hk_count: &HKCount,
    hyperkmers: &[String],
    minimizer: &str,
    extended_hyperkmer_left: &str,
) -> Option<usize> {
    let (_change_orientation, canonical_extended_hyperkmer_left) =
        get_canonical_kmer(extended_hyperkmer_left);
    for (candidate_left_extended_hk, _candidate_right_extended_hk, _count) in
        hk_count.get_iter(minimizer)
    {
        if canonical_extended_hyperkmer_left == hyperkmers[candidate_left_extended_hk.0] {
            return Some(candidate_left_extended_hk.0);
        }
    }
    None
}

fn get_start_of_minimizer_in_superkmer(sk: &SuperKmerInfos) -> usize {
    let m = sk.minimizer.len();
    if sk.was_read_canonical {
        sk.start_of_minimizer_as_read - sk.start_of_superkmer_as_read
    } else {
        sk.superkmer.len() + sk.start_of_superkmer_as_read - sk.start_of_minimizer_as_read - m
    }
}

/// Searches for `extended_hyperkmer_right` in `hk_count[minimizer]`
/// minimizer is assumed to be in canonical form
fn get_extended_hyperkmer_right_id(
    hk_count: &HKCount,
    hyperkmers: &[String],
    minimizer: &str,
    extended_hyperkmer_right: &str,
) -> Option<usize> {
    let (_change_orientation, canonical_extended_hyperkmer_right) =
        get_canonical_kmer(extended_hyperkmer_right);
    for (_candidate_left_extended_hk, candidate_right_extended_hk, _count) in
        hk_count.get_iter(minimizer)
    {
        if canonical_extended_hyperkmer_right == hyperkmers[candidate_right_extended_hk.0] {
            return Some(candidate_right_extended_hk.0);
        }
    }
    None
}

/// Adds `new_hyperkmer` in `hyperkmers` and return its index
/// `new_hyperkmer` does not have to be in canonical form
fn add_new_hyperkmer(hyperkmers: &mut Vec<String>, new_hyperkmer: &str) -> usize {
    let id_hyperkmer = hyperkmers.len();
    let new_hyperkmer = get_canonical_kmer(new_hyperkmer).1;
    hyperkmers.push(new_hyperkmer);
    id_hyperkmer
}

/// Add a new entry in the hyperkmer counts.
/// Minimizer should already be in canonical form.
fn insert_new_entry_in_hyperkmer_count(
    hk_count: &mut HKCount,
    minimizer: &str,
    id_left_hk: (usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
    id_right_hk: (usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
    count: Count,
) {
    hk_count.insert(String::from(minimizer), (id_left_hk, id_right_hk, count));
}

// OPTIMIZE use bytes instead of chars ?
fn common_suffix_length(x: &str, y: &str) -> usize {
    let mut x_chars = x.chars().rev();
    let mut y_chars = y.chars().rev();
    let mut length = 0;

    while let (Some(xc), Some(yc)) = (x_chars.next(), y_chars.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }

    length
}

fn common_prefix_length(x: &str, y: &str) -> usize {
    let mut x_chars = x.chars();
    let mut y_chars = y.chars();
    let mut length = 0;

    while let (Some(xc), Some(yc)) = (x_chars.next(), y_chars.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }

    length
}

// get the reverse complement of a sequence depending on the boolean parameter
// genereic implementation over the revers complement function to make tests easier
fn get_rc_if_change_orientation_internal<ReversComplementFunction>(
    revcompfunc: ReversComplementFunction,
    seq: &str,
    change_orientation: bool,
) -> String
where
    ReversComplementFunction: Fn(&str) -> String,
{
    if change_orientation {
        revcompfunc(seq)
    } else {
        String::from(seq)
    }
}

/// get the reverse complement of a sequence depending on the boolean parameter
fn get_rc_if_change_orientation(seq: &str, change_orientation: bool) -> String {
    get_rc_if_change_orientation_internal(reverse_complement, seq, change_orientation)
}

fn search_exact_hyperkmer_match(
    hyperkmers: &[String],
    left_hk: &str,
    right_hk: &str,
    candidate_left_ext_hk_metadata: &(usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
    candidate_right_ext_hk_metadata: &(usize, usize, PosOfHyperkmerInExtHyperkmer, bool),
) -> bool {
    // get sequences as they would appear if the current superkmer was canonical
    let candidate_left_ext_hk = get_rc_if_change_orientation(
        &hyperkmers[candidate_left_ext_hk_metadata.0],
        candidate_left_ext_hk_metadata.3,
    );
    let candidate_right_ext_hk = get_rc_if_change_orientation(
        &hyperkmers[candidate_right_ext_hk_metadata.0],
        candidate_right_ext_hk_metadata.3,
    );

    // extract candidate hyperkmers
    let candidate_left_hyperkmer = {
        match candidate_left_ext_hk_metadata.2 {
            PosOfHyperkmerInExtHyperkmer::Start => {
                &candidate_left_ext_hk[0..candidate_left_ext_hk_metadata.1]
            }
            PosOfHyperkmerInExtHyperkmer::End => {
                let end = candidate_left_ext_hk.len();
                &candidate_left_ext_hk[end - candidate_left_ext_hk_metadata.1..end]
            }
        }
    };
    let candidate_right_hyperkmer = {
        match candidate_right_ext_hk_metadata.2 {
            PosOfHyperkmerInExtHyperkmer::Start => {
                &candidate_right_ext_hk[0..candidate_right_ext_hk_metadata.1]
            }
            PosOfHyperkmerInExtHyperkmer::End => {
                let end = candidate_right_ext_hk.len();
                &candidate_right_ext_hk[end - candidate_right_ext_hk_metadata.1..end]
            }
        }
    };

    let match_left = candidate_left_hyperkmer == left_hk;
    let match_right = candidate_right_hyperkmer == right_hk;

    match_left && match_right
}

// TODO find a better name for the first stage function
fn first_stage(
    sequences: &Vec<&str>,
    k: usize,
    m: usize,
    threshold: Count,
) -> (SuperKmerCounts, HKCount, Vec<String>) {
    let mut sk_count = SuperKmerCounts::new();
    let mut hk_count: HKCount = MashMap::new();
    let mut hyperkmers: Vec<String> = Vec::new();

    for sequence in sequences {
        let start_superkmers = Instant::now();
        let superkmers = compute_superkmers_linear(sequence, k, m);
        println!(
            "super kmers computed in {} milliseconds",
            start_superkmers.elapsed().as_millis()
        );

        for triplet in superkmers.windows(3) {
            let (previous_sk, current_sk, next_sk) = (&triplet[0], &triplet[1], &triplet[2]);

            let (previous_sk, next_sk) = if current_sk.was_read_canonical {
                (previous_sk, next_sk)
            } else {
                (next_sk, previous_sk)
            };
            // now, current_sk.superkmer[0] is close to the left neighbour

            // compute now if the previous and/or next superkmer is solid
            // (we do it here and now because the incoming instruction `increase_count_superkmer(current_sk)`
            // can increase their count as well if they are the same superkmer)
            let previous_sk_is_solid = sk_count.get_count_superkmer(previous_sk) >= threshold;
            let next_sk_is_solid = sk_count.get_count_superkmer(next_sk) >= threshold;

            // TODO stocker les counts pour ne pas les recalculer
            let current_count = sk_count.increase_count_superkmer(current_sk);

            // chain of comparisons ahead, but I don't need to be exhaustive and I find it good as it is
            // so I tell clippy to shup up
            #[allow(clippy::comparison_chain)]
            if current_count == threshold {
                let (left_extended_hk, right_extended_hk) =
                    get_left_and_rigth_extended_hk(previous_sk, current_sk, next_sk);
                let (left_hk, right_hk) = get_left_and_rigth_hk(previous_sk, current_sk, next_sk);
                assert!(left_extended_hk.0.contains(&left_hk));
                assert!(left_extended_hk.0.len() < k);
                assert!(right_extended_hk.0.len() < k);
                assert!(right_extended_hk.0.contains(&right_hk));
                let canonical_extended_left_hk = get_canonical_kmer(&left_extended_hk.0);
                let canonical_extended_right_hk = get_canonical_kmer(&right_extended_hk.0);

                // left_hk and right_hk are the left and right hyperkmer
                // as we would see them if the minimizer was in canonical form in the read

                // OPTIMIZE maybe it is posssible to call get_hyperkmer_{left, right}_id and ignore get_count_superkmer
                // OPTIMIZE of even better: access the count of sk in streaming, so that no recomputation is needed
                // println!("{}", sk_count.get_count_superkmer(previous_sk));
                let id_left_hk = if previous_sk_is_solid {
                    if previous_sk.was_read_canonical == current_sk.was_read_canonical {
                        get_extended_hyperkmer_right_id(
                            &hk_count,
                            &hyperkmers,
                            &previous_sk.minimizer,
                            &canonical_extended_left_hk.1,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        get_extended_hyperkmer_left_id(
                            &hk_count,
                            &hyperkmers,
                            &previous_sk.minimizer,
                            &canonical_extended_left_hk.1,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &left_extended_hk.0)
                };
                let id_right_hk = if next_sk_is_solid {
                    if current_sk.was_read_canonical == next_sk.was_read_canonical {
                        get_extended_hyperkmer_left_id(
                            &hk_count,
                            &hyperkmers,
                            &next_sk.minimizer,
                            &right_extended_hk.0,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        get_extended_hyperkmer_right_id(
                            &hk_count,
                            &hyperkmers,
                            &next_sk.minimizer,
                            &right_extended_hk.0,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &right_extended_hk.0)
                };

                // we have two ids (left and rigth) of extended hyperkmers containing our left and right hyperkemr
                // let's get their orientation wrt to the orientation of the minimizer

                let left_change_orientation = !canonical_extended_left_hk.0;
                let right_change_orientation = !canonical_extended_right_hk.0;

                let candidate_left_ext_hk = if left_change_orientation {
                    reverse_complement(&hyperkmers[id_left_hk])
                } else {
                    hyperkmers[id_left_hk].clone()
                };
                let candidate_right_ext_hk = if right_change_orientation {
                    reverse_complement(&hyperkmers[id_right_hk])
                } else {
                    hyperkmers[id_right_hk].clone()
                };

                assert!(candidate_left_ext_hk == left_extended_hk.0);
                assert!(candidate_right_ext_hk == right_extended_hk.0);
                assert!(hyperkmers[id_left_hk] == get_canonical_kmer(&left_extended_hk.0).1);
                assert!(hyperkmers[id_right_hk] == get_canonical_kmer(&right_extended_hk.0).1);

                insert_new_entry_in_hyperkmer_count(
                    &mut hk_count,
                    &current_sk.minimizer,
                    (
                        id_left_hk,
                        left_extended_hk.2,
                        left_extended_hk.1,
                        left_change_orientation,
                    ),
                    (
                        id_right_hk,
                        right_extended_hk.2,
                        right_extended_hk.1,
                        right_change_orientation,
                    ),
                    current_count,
                );
            } else if current_count > threshold {
                let (left_hk, right_hk) = get_left_and_rigth_hk(previous_sk, current_sk, next_sk);
                let mut found = false;
                for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, count_hk) in
                    hk_count.get_mut_iter(&current_sk.minimizer)
                {
                    let is_exact_match = search_exact_hyperkmer_match(
                        &hyperkmers,
                        &left_hk,
                        &right_hk,
                        candidate_left_ext_hk_metadata,
                        candidate_right_ext_hk_metadata,
                    );
                    if is_exact_match {
                        *count_hk += 1;
                        found = true;
                        break;
                    }
                }
                if !found {
                    let mut new_left_and_right_metadata = None;
                    for (
                        candidate_left_ext_hk_metadata,
                        candidate_right_ext_hk_metadata,
                        _count_hk,
                    ) in hk_count.get_mut_iter(&current_sk.minimizer)
                    {
                        // get sequences as they would appear if the current superkmer was canonical
                        let candidate_left_ext_hk = if candidate_left_ext_hk_metadata.3 {
                            reverse_complement(&hyperkmers[candidate_left_ext_hk_metadata.0])
                        } else {
                            hyperkmers[candidate_left_ext_hk_metadata.0].clone()
                        };
                        let candidate_right_ext_hk = if candidate_right_ext_hk_metadata.3 {
                            reverse_complement(&hyperkmers[candidate_right_ext_hk_metadata.0])
                        } else {
                            hyperkmers[candidate_right_ext_hk_metadata.0].clone()
                        };

                        let match_start_left = candidate_left_ext_hk.starts_with(&left_hk);
                        let match_end_left = candidate_left_ext_hk.ends_with(&left_hk);

                        let match_start_right = candidate_right_ext_hk.starts_with(&right_hk);
                        let match_end_right = candidate_right_ext_hk.ends_with(&right_hk);

                        // TODO code duplication
                        if match_start_left && match_start_right {
                            let metadata_to_insert_left = (
                                candidate_left_ext_hk_metadata.0,
                                left_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::Start,
                                candidate_left_ext_hk_metadata.3,
                            );
                            let metadata_to_insert_right = (
                                candidate_right_ext_hk_metadata.0,
                                right_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::Start,
                                candidate_right_ext_hk_metadata.3,
                            );
                            new_left_and_right_metadata =
                                Some((metadata_to_insert_left, metadata_to_insert_right));
                            break;
                        } else if match_end_left && match_start_right {
                            let metadata_to_insert_left = (
                                candidate_left_ext_hk_metadata.0,
                                left_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::End,
                                candidate_left_ext_hk_metadata.3,
                            );
                            let metadata_to_insert_right = (
                                candidate_right_ext_hk_metadata.0,
                                right_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::Start,
                                candidate_right_ext_hk_metadata.3,
                            );
                            new_left_and_right_metadata =
                                Some((metadata_to_insert_left, metadata_to_insert_right));
                            break;
                        } else if match_start_left && match_end_right {
                            let metadata_to_insert_left = (
                                candidate_left_ext_hk_metadata.0,
                                left_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::Start,
                                candidate_left_ext_hk_metadata.3,
                            );
                            let metadata_to_insert_right = (
                                candidate_right_ext_hk_metadata.0,
                                right_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::End,
                                candidate_right_ext_hk_metadata.3,
                            );
                            new_left_and_right_metadata =
                                Some((metadata_to_insert_left, metadata_to_insert_right));
                            break;
                        } else if match_end_left && match_end_right {
                            let metadata_to_insert_left = (
                                candidate_left_ext_hk_metadata.0,
                                left_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::End,
                                candidate_left_ext_hk_metadata.3,
                            );
                            let metadata_to_insert_right = (
                                candidate_right_ext_hk_metadata.0,
                                right_hk.len(),
                                PosOfHyperkmerInExtHyperkmer::End,
                                candidate_right_ext_hk_metadata.3,
                            );
                            new_left_and_right_metadata =
                                Some((metadata_to_insert_left, metadata_to_insert_right));
                            break;
                        }
                    }
                    // TODO error message here
                    let (metadata_to_insert_left, metadata_to_insert_right) =
                        new_left_and_right_metadata.unwrap();
                    insert_new_entry_in_hyperkmer_count(
                        &mut hk_count,
                        &current_sk.minimizer,
                        metadata_to_insert_left,
                        metadata_to_insert_right,
                        current_count,
                    );
                    // ensure the hyperkmer was correctly inserted
                    assert!(search_exact_hyperkmer_match(
                        &hyperkmers,
                        &left_hk,
                        &right_hk,
                        &metadata_to_insert_left,
                        &metadata_to_insert_right
                    ));
                }
            }
        }
    }
    (sk_count, hk_count, hyperkmers)
}

// TODO find a better name for the second stage function
fn second_stage(
    sk_count: &mut SuperKmerCounts,
    hk_count: &mut HKCount,
    hyperkmers: &[String],
    sequences: &Vec<&str>,
    k: usize,
    m: usize,
    threshold: Count,
) -> (TruncatedHKCount, HashMap<Minimizer, Count>) {
    let mut truncated_hk = TruncatedHKCount::new();
    // let mut truncated_hk: MashMap<Minimizer, (usize, usize)> = MashMap::new();
    let mut discarded_minimizers: HashMap<Minimizer, Count> = HashMap::new();

    for sequence in sequences {
        let superkmers = compute_superkmers_linear(sequence, k, m);
        for superkmer in &superkmers {
            if sk_count.get_count_superkmer(superkmer) >= threshold {
                continue;
            }

            let minimizer = &superkmer.minimizer;

            if !hk_count.contains_key(minimizer) {
                // increase count, set to 1 if it was 0
                let count = discarded_minimizers
                    .entry(minimizer.clone())
                    .and_modify(|counter| {
                        *counter = counter.saturating_add(1);
                    })
                    .or_insert(1);
                if *count == threshold {
                    warn!(
                        "minimizer {} of superkmer {} is found {} times but its hyperkmer is not",
                        minimizer, superkmer.superkmer, count
                    );
                }
            }

            let start_of_minimizer_in_sk = get_start_of_minimizer_in_superkmer(superkmer);
            let left_sk = &superkmer.superkmer[0..(start_of_minimizer_in_sk + m - 1)];
            let right_sk =
                &superkmer.superkmer[(start_of_minimizer_in_sk + 1)..superkmer.superkmer.len()];

            let mut match_size = k - 1;
            let mut match_metadata = None;

            // Search for the maximal preffix/suffix in hyperkmers.
            // Because we only have our left/right *partial* hyperkmer,
            // we do not know if they are in the same orientation as their "complete" canonical form.
            // Let's search for both orientation then.
            for (candidate_left_ext_hk_metadata, candidate_right_ext_hk_metadata, _count) in
                hk_count.get_iter(minimizer)
            {
                let candidate_left_ext_hk = get_rc_if_change_orientation(
                    &hyperkmers[candidate_left_ext_hk_metadata.0],
                    candidate_left_ext_hk_metadata.3,
                );
                let candidate_right_ext_hk = get_rc_if_change_orientation(
                    &hyperkmers[candidate_right_ext_hk_metadata.0],
                    candidate_right_ext_hk_metadata.3,
                );

                // extract candidate hyperkmers
                let candidate_left_hyperkmer = {
                    match candidate_left_ext_hk_metadata.2 {
                        PosOfHyperkmerInExtHyperkmer::Start => {
                            &candidate_left_ext_hk[0..candidate_left_ext_hk_metadata.1]
                        }
                        PosOfHyperkmerInExtHyperkmer::End => {
                            let end = candidate_left_ext_hk.len();
                            &candidate_left_ext_hk[end - candidate_left_ext_hk_metadata.1..end]
                        }
                    }
                };
                let candidate_right_hyperkmer = {
                    match candidate_right_ext_hk_metadata.2 {
                        PosOfHyperkmerInExtHyperkmer::Start => {
                            &candidate_right_ext_hk[0..candidate_right_ext_hk_metadata.1]
                        }
                        PosOfHyperkmerInExtHyperkmer::End => {
                            let end = candidate_right_ext_hk.len();
                            &candidate_right_ext_hk[end - candidate_right_ext_hk_metadata.1..end]
                        }
                    }
                };

                let len_current_match_left =
                    common_suffix_length(left_sk, candidate_left_hyperkmer);
                let len_current_match_right =
                    common_prefix_length(right_sk, candidate_right_hyperkmer);
                let current_match_size = len_current_match_left + len_current_match_right;

                if current_match_size > match_size {
                    match_size = current_match_size;
                    // extract candidate hyperkmers
                    let left_end = {
                        match candidate_left_ext_hk_metadata.2 {
                            PosOfHyperkmerInExtHyperkmer::Start => candidate_left_ext_hk_metadata.1,
                            PosOfHyperkmerInExtHyperkmer::End => candidate_left_ext_hk.len(),
                        }
                    };
                    let left_start = left_end - len_current_match_left;

                    let right_start = {
                        match candidate_right_ext_hk_metadata.2 {
                            PosOfHyperkmerInExtHyperkmer::Start => 0,
                            PosOfHyperkmerInExtHyperkmer::End => {
                                let end = candidate_right_ext_hk.len();
                                end - candidate_right_ext_hk_metadata.1
                            }
                        }
                    };
                    let right_end = right_start + len_current_match_right;

                    match_metadata = Some((
                        (
                            candidate_left_ext_hk_metadata.0,
                            left_start,
                            left_end,
                            candidate_left_ext_hk_metadata.3,
                        ),
                        (
                            candidate_right_ext_hk_metadata.0,
                            right_start,
                            right_end,
                            candidate_right_ext_hk_metadata.3,
                        ),
                    ));
                }
            }

            if let Some(metadata) = match_metadata {
                truncated_hk.insert(minimizer.clone(), (metadata.0, metadata.1, 1));
            }
        }
    }
    (truncated_hk, discarded_minimizers)
}

fn index_hyperkmers(
    k: usize,
    m: usize,
    threshold: Count,
    sequences: &Vec<&str>,
) -> (
    SuperKmerCounts,
    HKCount,
    Vec<String>,
    TruncatedHKCount,
    HashMap<String, u16>,
) {
    let start_fisrt_step = Instant::now();
    let (mut sk_count, mut hk_count, hyperkmers) = first_stage(sequences, k, m, threshold);
    println!(
        "time first stage: {} milliseconds",
        start_fisrt_step.elapsed().as_millis()
    );
    let start_second_stage = Instant::now();
    let (truncated_hk, discarded_minimizers) = second_stage(
        &mut sk_count,
        &mut hk_count,
        &hyperkmers,
        sequences,
        k,
        m,
        threshold,
    );
    println!(
        "time second stage: {} milliseconds",
        start_second_stage.elapsed().as_millis()
    );
    (
        sk_count,
        hk_count,
        hyperkmers,
        truncated_hk,
        discarded_minimizers,
    )
}

fn main() {
    let sequences: Vec<String> = fasta_reads("data/U00096.3.fasta")
        .unwrap()
        .map(|rcstring| rcstring.to_string())
        .collect();
    let sequences: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();

    let k = 31;
    let m = 20;
    let threshold = 2;

    let (_sk_count, hk_count, hyperkmers, _truncated_hk, _discarded_minimizers) =
        index_hyperkmers(k, m, threshold, &sequences);

    let kmer_test = "CGCGAGGAGCTGGCCGAGGTGGATGTGGACTGGCTGATCGCCGAGCGCCCCGGCAAGGTAAGAACCTTGAAACAGCATCCACGCAAGAACAAAACGGCCA";

    // search(&hk_count, &hyperkmers, kmer_test, k, m);
    // println!("{:?}", sk);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suffix_preffix() {
        let x = "aazerty______";
        let y = "aazerty-___---------";
        assert_eq!(common_prefix_length(x, y), 7);
        assert_eq!(common_prefix_length(y, x), 7);

        assert_eq!(common_prefix_length("x", "y"), 0);
        assert_eq!(common_prefix_length("", ""), 0);
        assert_eq!(common_prefix_length("xyyyyyyyyyy", "xyy"), 3);

        // suffix
        let x = "aazerty_____-erc---------";
        let y = "aazerty-___erc---------";
        assert_eq!(common_suffix_length(x, y), 12);
        assert_eq!(common_suffix_length(y, x), 12);

        assert_eq!(common_suffix_length("x", "y"), 0);
        assert_eq!(common_suffix_length("", ""), 0);
        assert_eq!(common_suffix_length("xyyyyyyyyyy", "xyy"), 2);
    }

    #[test]
    fn test_get_left_and_rigth_extended_hk() {
        // in the sequence CAGATGGTTCAACCCTTAAGTTAGCGCTTATGGGATCAC
        //     previous sk :              CCTTAAGTTAGCGCTTATGGGATCAC (GTTAGCGCTT)
        //     current sk :            ACCCTTAAGTTAGCGCTTATG (CCCTTAAGTT)
        //        next sk : CAGATGGTTCAACCCTTAAGTTAGCGCTTA (AACCCTTAAG)

        //    SuperKmerInfos { superkmer: "GTGATCCCATAAGCGCTAACTTAAGG", minimizer: "AAGCGCTAAC", was_read_canonical: false, start_of_minimizer_as_read: 16715, start_of_superkmer_as_read: 16709 }
        //    SuperKmerInfos { superkmer: "CATAAGCGCTAACTTAAGGGT", minimizer: "AACTTAAGGG", was_read_canonical: false, start_of_minimizer_as_read: 16708, start_of_superkmer_as_read: 16707 }
        //    SuperKmerInfos { superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA", minimizer: "AACCCTTAAG", was_read_canonical: true, start_of_minimizer_as_read: 16706, start_of_superkmer_as_read: 16696 }
        //    ("ACTTAAGGGT", "TAAGCGCTAACTTAAGGGT")
        //    left extended hk: ("AGCGCTAACTTAAGG", 15)
        //    right extended hk: ("TAAGCGCTAACTTAAGGGT", 10)
        //    TAAGCGCTAACTTAAGGGTTGAACCATCTG
        let previous_sk = SuperKmerInfos {
            superkmer: "GTGATCCCATAAGCGCTAACTTAAGG".into(),
            minimizer: "AAGCGCTAAC".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 16715,
            start_of_superkmer_as_read: 16709,
        };
        let current_sk = SuperKmerInfos {
            superkmer: "CATAAGCGCTAACTTAAGGGT".into(),
            minimizer: "AACTTAAGGG".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 16708,
            start_of_superkmer_as_read: 16707,
        };
        let next_sk = SuperKmerInfos {
            superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA".into(),
            minimizer: "AACCCTTAAG".into(),
            was_read_canonical: true,
            start_of_minimizer_as_read: 16706,
            start_of_superkmer_as_read: 16696,
        };
        let (l, r) = get_left_and_rigth_extended_hk(&previous_sk, &current_sk, &next_sk);
        //   left: "TAAGCGCTAACTTAAGGGT"
        //  right: "ACCCTTAAGTTAGCGCTTA"
        assert_eq!(r.0, "TAAGCGCTAACTTAAGGGT");
        assert_eq!(l.0, "CATAAGCGCTAACTTAAGG");
        //    sk right: ACCCTTAAGTTAGCGCTTA / TAAGCGCTAACTTAAGGGT
        //    sk left: CCTTAAGTTAGCGCTTATG / CATAAGCGCTAACTTAAGG
    }

    #[test]
    fn test_get_hyperkmer_left_and_rigth_id() {
        use rand::{distributions::Alphanumeric, Rng}; // 0.8
        let mut hk_count: HKCount = MashMap::new();
        let mut hyperkmers: Vec<String> = Vec::new();

        let nb_insertions = 10000;

        // random insertions in hyperkmers
        for _ in 0..nb_insertions {
            let s: String = rand::thread_rng()
                .sample_iter(&Alphanumeric)
                .take(100)
                .map(char::from)
                .collect();
            hyperkmers.push(get_canonical_kmer(&s).1);
        }

        // random minimizers
        let mut minimizers = Vec::new();
        for _ in 0..(hyperkmers.len() - 1) {
            let minimizer: String = rand::thread_rng()
                .sample_iter(&Alphanumeric)
                .take(20)
                .map(char::from)
                .collect();
            minimizers.push(get_canonical_kmer(&minimizer).1);
        }

        // link minimizer and hyperkmers
        // start by inserting random values
        let mut rng = rand::thread_rng();
        for _ in 0..3 {
            for minimizer in &minimizers {
                let random_left = rng.gen_range(0..hyperkmers.len());
                let random_rigth = rng.gen_range(0..hyperkmers.len());
                let random_left_overlap = rng.gen_range(0..hyperkmers.len());
                let random_rigth_overlap = rng.gen_range(0..hyperkmers.len());
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(
                    minimizer.clone(),
                    (
                        (
                            random_left,
                            random_left_overlap,
                            PosOfHyperkmerInExtHyperkmer::Start,
                            true,
                        ),
                        (
                            random_rigth,
                            random_rigth_overlap,
                            PosOfHyperkmerInExtHyperkmer::Start,
                            true,
                        ),
                        random_count,
                    ),
                );
            }
        }

        // then link minimizer and hyperkmers
        for (i, minimizer) in minimizers.iter().enumerate() {
            let random_count: u16 = rng.gen_range(0..hyperkmers.len()) as u16;
            let overlap = minimizer.len();
            hk_count.insert(
                minimizer.clone(),
                (
                    (i, overlap, PosOfHyperkmerInExtHyperkmer::Start, true),
                    (i + 1, overlap, PosOfHyperkmerInExtHyperkmer::Start, true),
                    random_count,
                ),
            );
        }

        // add another random values
        for _ in 0..10 {
            for minimizer in &minimizers {
                let random_left = rng.gen_range(0..hyperkmers.len());
                let random_rigth = rng.gen_range(0..hyperkmers.len());
                let random_left_overlap = rng.gen_range(0..hyperkmers.len());
                let random_rigth_overlap = rng.gen_range(0..hyperkmers.len());
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(
                    minimizer.clone(),
                    (
                        (
                            random_left,
                            random_left_overlap,
                            PosOfHyperkmerInExtHyperkmer::Start,
                            true,
                        ),
                        (
                            random_rigth,
                            random_rigth_overlap,
                            PosOfHyperkmerInExtHyperkmer::Start,
                            true,
                        ),
                        random_count,
                    ),
                );
            }
        }

        // among all the random values, we are still able to get back our real data
        let hk_id =
            get_extended_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizers[5], &hyperkmers[5]);
        assert_eq!(hk_id, Some(5));
        let hk_id =
            get_extended_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizers[5], &hyperkmers[6]);
        assert_eq!(hk_id, Some(6));

        // query a non existant minimiser
        let minimizer: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(20)
            .map(char::from)
            .collect();
        let hk_id =
            get_extended_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[5]);
        assert_eq!(hk_id, None);
        let hk_id =
            get_extended_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[6]);
        assert_eq!(hk_id, None);

        // query an existing minimizer, but with a random hyperkmer
        let hyperkmer: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(100)
            .map(char::from)
            .collect();
        let hk_id = get_extended_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmer);
        assert_eq!(hk_id, None);
        let hk_id = get_extended_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmer);
        assert_eq!(hk_id, None);

        // query an existing minimizer, but with a random existing hyperkmer (TODO probability of failure is small but not null)
        let hk_id =
            get_extended_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[3]);
        assert_eq!(hk_id, None);
        let hk_id =
            get_extended_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[4]);
        assert_eq!(hk_id, None);
    }

    // #[test]
    // fn test_suf_of_x_is_pref_of_y() {
    //     assert_eq!(suf_of_x_is_pref_of_y("abcdefghij", "fghijyhy"), 5);
    // }
}
