use fastxgz::fasta_reads;
use log::warn;
use mashmap::MashMap;
use std::collections::HashMap;
use std::time::Instant;

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;
// type SKCount = MashMap<Minimizer, (HashSuperKmer, Count)>;
// canonical minimmizer -> (id of left hk, orientation flag), ((index of right hyperkmer, orientation flag), count)
// orientation flag: true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer
type HKCount = MashMap<Minimizer, ((usize, usize, bool), (usize, usize, bool), Count)>;

mod brrr_minimizers;
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
    println!(
        "get_extended_hyperkmer_left_id, searching for {}",
        canonical_extended_hyperkmer_left
    );
    for (candidate_left_extended_hk, _candidate_right_extended_hk, _count) in
        hk_count.get_iter(minimizer)
    {
        println!("candidate: {:?}", hyperkmers[candidate_left_extended_hk.0]);
        if canonical_extended_hyperkmer_left == hyperkmers[candidate_left_extended_hk.0] {
            return Some(candidate_left_extended_hk.0);
        }
    }
    None
}

fn get_left_part_of_canonical_superkmer(sk: &SuperKmerInfos) {}
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
    println!(
        "get_extended_hyperkmer_right_id, searching for {}",
        canonical_extended_hyperkmer_right
    );
    if !hk_count.contains_key(minimizer) {
        println!("table does not have key {}", minimizer);
    }
    for (_candidate_left_extended_hk, candidate_right_extended_hk, _count) in
        hk_count.get_iter(minimizer)
    {
        println!(
            "searching for ext hk {}, current ext hk is{}",
            canonical_extended_hyperkmer_right, hyperkmers[candidate_right_extended_hk.0]
        );
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
    id_left_hk: (usize, usize, bool),
    id_right_hk: (usize, usize, bool),
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

/// Would return 6 for ACTGGC**TGCGTAGC** **TGCGTAGC**TGCA
/// Would return 6
fn suf_of_x_is_pref_of_y(x: &str, y: &str) -> usize {
    let max_overlap = x.len().min(y.len());

    for i in (0..=max_overlap).rev() {
        if x.ends_with(&y[..i]) {
            return i;
        }
    }

    0
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

fn get_left_and_rigth_hk(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> (String, String) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

    let (start_of_minimizer_in_sk, distance_to_left, distance_to_right) =
        if current_sk.was_read_canonical {
            assert!(current_sk.start_of_minimizer_as_read > previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read < next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk =
                current_sk.start_of_minimizer_as_read - current_sk.start_of_superkmer_as_read;
            let distance_to_left =
                current_sk.start_of_minimizer_as_read - previous_sk.start_of_minimizer_as_read;
            let distance_to_right =
                next_sk.start_of_minimizer_as_read - current_sk.start_of_minimizer_as_read;
            (
                start_of_minimizer_in_sk,
                distance_to_left,
                distance_to_right,
            )
        } else {
            // we work on the rc of the superkmer
            // as distances are given as if the superkmer is canonical in the read, we must do some math
            assert!(current_sk.start_of_minimizer_as_read < previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read > next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk = current_sk.superkmer.len()
                - current_sk.start_of_minimizer_as_read
                + current_sk.start_of_superkmer_as_read
                - m;
            let distance_to_left =
                previous_sk.start_of_minimizer_as_read - current_sk.start_of_minimizer_as_read;
            let distance_to_right =
                current_sk.start_of_minimizer_as_read - next_sk.start_of_minimizer_as_read;
            (
                start_of_minimizer_in_sk,
                distance_to_left,
                distance_to_right,
            )
        };

    let start_of_left_hk = start_of_minimizer_in_sk - distance_to_left + 1;
    let end_of_left_hk = start_of_minimizer_in_sk + m - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;
    let end_of_right_hk = start_of_minimizer_in_sk + distance_to_right + m - 1;

    // this is the sequence in between minimizers, including m-1 bases of the minimizers
    let left_hk = &current_sk.superkmer[start_of_left_hk..end_of_left_hk];
    let right_hk = &current_sk.superkmer[start_of_right_hk..end_of_right_hk];

    (String::from(left_hk), String::from(right_hk))
}

/// return the left and right extended hyperkmers of the currrent superkmer
/// extended hyperkmers are, in fact, part of superkmers with the length of inclusion of the hyperkmer in it
/// returned sequences are not in canonical form, but are sequences as they would appear in the read
/// if the current superkmer was in its canonical form in the read
/// returned tuple is (left, right)
/// left is (seq, i) st hyperkmer = seq[i:]
/// right is (seq, i) st hyperkmer = seq[:-i]
fn get_left_and_rigth_extended_hk(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> ((String, usize), (String, usize)) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

    // get the starting position of the minimizer in the superkmer
    let get_start_of_minimizer = |sk: &SuperKmerInfos| {
        if sk.was_read_canonical {
            sk.start_of_minimizer_as_read - sk.start_of_superkmer_as_read
        } else {
            sk.superkmer.len() + sk.start_of_superkmer_as_read - sk.start_of_minimizer_as_read - m
        }
    };

    let (start_of_minimizer_in_sk, distance_to_left, distance_to_right) =
        if current_sk.was_read_canonical {
            assert!(current_sk.start_of_minimizer_as_read > previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read < next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk =
                current_sk.start_of_minimizer_as_read - current_sk.start_of_superkmer_as_read;
            let distance_to_left =
                current_sk.start_of_minimizer_as_read - previous_sk.start_of_minimizer_as_read;
            let distance_to_right =
                next_sk.start_of_minimizer_as_read - current_sk.start_of_minimizer_as_read;
            (
                start_of_minimizer_in_sk,
                distance_to_left,
                distance_to_right,
            )
        } else {
            // we work on the rc of the superkmer
            // as distances are given as if the superkmer is canonical in the read, we must do some math
            assert!(current_sk.start_of_minimizer_as_read < previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read > next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk = current_sk.superkmer.len()
                + current_sk.start_of_superkmer_as_read
                - current_sk.start_of_minimizer_as_read
                - m;
            let distance_to_left =
                previous_sk.start_of_minimizer_as_read - current_sk.start_of_minimizer_as_read;
            let distance_to_right =
                current_sk.start_of_minimizer_as_read - next_sk.start_of_minimizer_as_read;
            (
                start_of_minimizer_in_sk,
                distance_to_left,
                distance_to_right,
            )
        };

    let start_of_left_hk = start_of_minimizer_in_sk - distance_to_left + 1;
    let end_of_left_hk = start_of_minimizer_in_sk + m - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;
    let end_of_right_hk = start_of_minimizer_in_sk + distance_to_right + m - 1;

    // this is the sequence in between minimizers, including m-1 bases of the minimizers
    let left_hk = &current_sk.superkmer[start_of_left_hk..end_of_left_hk];
    let right_hk = &current_sk.superkmer[start_of_right_hk..end_of_right_hk];
    // this is the left and rigth part of the superkmer
    let current_left_sk = &current_sk.superkmer[0..end_of_left_hk];
    let current_right_sk = &current_sk.superkmer[start_of_right_hk..current_sk.superkmer.len()];

    let previous_right_sk = if previous_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(
            &previous_sk.superkmer
                [get_start_of_minimizer(previous_sk) + 1..previous_sk.superkmer.len()],
        )
    } else {
        reverse_complement(&previous_sk.superkmer[0..get_start_of_minimizer(previous_sk) + m - 1])
    };
    let next_left_sk = if next_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(&next_sk.superkmer[0..get_start_of_minimizer(next_sk) + m - 1])
    } else {
        reverse_complement(
            &next_sk.superkmer[get_start_of_minimizer(next_sk) + 1..next_sk.superkmer.len()],
        )
    };
    // TODO remove checks
    let extended_left_sk = if current_left_sk.len() > previous_right_sk.len() {
        assert!(current_left_sk.contains(&previous_right_sk));
        assert!(current_left_sk.ends_with(&previous_right_sk));
        current_left_sk.into()
    } else {
        assert!(previous_right_sk.contains(current_left_sk));
        assert!(previous_right_sk.starts_with(current_left_sk));
        previous_right_sk.clone()
    };
    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        assert!(current_right_sk.contains(&next_left_sk));
        assert!(current_right_sk.starts_with(&next_left_sk));
        current_right_sk.into()
    } else {
        assert!(next_left_sk.contains(current_right_sk));
        assert!(next_left_sk.ends_with(current_right_sk));
        next_left_sk.clone()
    };
    return (
        (
            extended_left_sk,
            std::cmp::min(current_left_sk.len(), previous_right_sk.len()),
        ),
        (
            extended_right_sk,
            std::cmp::min(current_right_sk.len(), next_left_sk.len()),
        ),
    );
    // now, let's compute left and right hyperkmers and check that the hyperkmers are in them

    let overlap_right = common_suffix_length(right_hk, &next_left_sk);
    let larger_right_sk = if right_hk.len() > next_left_sk.len() {
        right_hk
    } else {
        &next_left_sk
    };

    let overlap_left = suf_of_x_is_pref_of_y(left_hk, &previous_right_sk);
    let larger_left_sk = if left_hk.len() > previous_right_sk.len() {
        left_hk
    } else {
        &previous_right_sk
    };

    println!("{:?}", overlap_right);
    println!("{:?}\n{:?}\n{:?}", previous_sk, current_sk, next_sk);
    println!("{:?}", (right_hk, &next_left_sk));
    // sanity checks
    // overlaps are longer or equal than len(minimizers)-1
    assert!(overlap_right >= current_sk.minimizer.len() - 1);
    assert!(overlap_left >= current_sk.minimizer.len() - 1);
    // one sequence should be completely covered by the overlap
    assert_eq!(
        overlap_right,
        std::cmp::min(right_hk.len(), next_left_sk.len())
    );
    assert_eq!(
        overlap_left,
        std::cmp::min(left_hk.len(), previous_right_sk.len())
    );
    // each overlap contains both sequences
    assert!(larger_right_sk.contains(right_hk));
    assert!(larger_right_sk.contains(&next_left_sk));
    assert!(larger_left_sk.contains(left_hk));
    assert!(larger_left_sk.contains(&previous_right_sk));

    (
        (String::from(larger_left_sk), overlap_left),
        (String::from(larger_right_sk), overlap_right),
    )
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

            // compute now if the previosu and/or next superkmer is solid
            // (we do it now because `increase_count_superkmer(current_sk)`
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
                    println!(
                        "sk_count contains_minimizer: {}",
                        sk_count.contains_minimizer(&previous_sk.minimizer)
                    );
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

                let left_change_orientation =
                    canonical_extended_left_hk.0 != current_sk.was_read_canonical;
                let right_change_orientation =
                    canonical_extended_right_hk.0 != current_sk.was_read_canonical;

                insert_new_entry_in_hyperkmer_count(
                    &mut hk_count,
                    &current_sk.minimizer,
                    (id_left_hk, left_extended_hk.1, left_change_orientation),
                    (id_right_hk, right_extended_hk.1, right_change_orientation),
                    current_count,
                );
            } else if current_count > threshold {
                let mut found = false;
                for (candidate_left_hk, candidate_right_hk, count_hk) in
                    hk_count.get_mut_iter(&current_sk.minimizer)
                {
                    let mut left_hk = hyperkmers[candidate_left_hk.0].clone();
                    if candidate_left_hk.2 {
                        left_hk = reverse_complement(&left_hk);
                    }

                    let mut right_hk = hyperkmers[candidate_right_hk.0].clone();
                    if candidate_right_hk.2 {
                        right_hk = reverse_complement(&right_hk);
                    }

                    if left_hk.ends_with(&left_hk) && right_hk.starts_with(&right_hk) {
                        *count_hk += 1;
                        found = true;
                        break;
                    }
                }
                if !found {
                    println!("error");
                    panic!();
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
) -> (
    MashMap<std::string::String, (usize, usize)>,
    HashMap<std::string::String, u16>,
) {
    let mut truncated_hk: MashMap<Minimizer, (usize, usize)> = MashMap::new();
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

            let start_of_minimizer_in_sk =
                superkmer.start_of_minimizer_as_read - superkmer.start_of_superkmer_as_read;
            let left_hk_of_sk = &superkmer.superkmer[0..(start_of_minimizer_in_sk + m - 1)];
            let right_hk_of_sk =
                &superkmer.superkmer[(start_of_minimizer_in_sk + 1)..superkmer.superkmer.len()];

            let mut match_left = k - 1;
            let mut match_right = k - 1;
            let mut id_match_left = None;
            let mut id_match_right = None;

            // Search for the maximal preffix/suffix in hyperkmers.
            // Because we only have our left/right *partial* hyperkmer,
            // we do not know if they are in the same orientation as their "complete" canonical form.
            // Let's search for both orientation then.
            for (id_left_hk, id_right_hk, _count) in hk_count.get_iter(minimizer) {
                let len_current_match_left = std::cmp::max(
                    common_suffix_length(left_hk_of_sk, &hyperkmers[id_left_hk.0]),
                    common_suffix_length(
                        &reverse_complement(left_hk_of_sk),
                        &hyperkmers[id_left_hk.0],
                    ),
                );

                if len_current_match_left > match_left {
                    match_left = len_current_match_left;
                    id_match_left = Some(id_left_hk);
                }

                let len_current_match_right = std::cmp::max(
                    common_prefix_length(right_hk_of_sk, &hyperkmers[id_right_hk.0]),
                    common_prefix_length(
                        &reverse_complement(right_hk_of_sk),
                        &hyperkmers[id_right_hk.0],
                    ),
                );
                if len_current_match_right > match_right {
                    match_right = len_current_match_right;
                    id_match_right = Some(id_right_hk);
                }
            }

            if let Some(id_match_left) = id_match_left {
                truncated_hk.insert(minimizer.clone(), (id_match_left.0, match_left));
            }

            if let Some(id_match_right) = id_match_right {
                truncated_hk.insert(minimizer.clone(), (id_match_right.0, match_left));
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
    MashMap<String, (usize, usize)>,
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
                        (random_left, random_left_overlap, true),
                        (random_rigth, random_rigth_overlap, true),
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
                ((i, overlap, true), (i + 1, overlap, true), random_count),
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
                        (random_left, random_left_overlap, true),
                        (random_rigth, random_rigth_overlap, true),
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

    #[test]
    fn test_suf_of_x_is_pref_of_y() {
        assert_eq!(suf_of_x_is_pref_of_y("abcdefghij", "fghijyhy"), 5);
    }
}
