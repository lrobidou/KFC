use log::warn;
use mashmap::MashMap;
use std::collections::HashMap;
use std::time::Instant;
use xxhash_rust::const_xxh3::xxh3_64;

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;
type SKCount = MashMap<Minimizer, (HashSuperKmer, Count)>;
// canonical minimmizer -> (id of left hk, orientation flag), ((index of right hyperkmer, orientation flag), count)
// orientation flag: true if the hyperkmer is in a DIFFERENT orientation that the canonical minimizer
type HKCount = MashMap<Minimizer, ((usize, bool), (usize, bool), Count)>;

mod brrr_minimizers;
mod superkmers;

use fastxgz::fasta_reads;
use superkmers::{
    compute_superkmers_linear, get_canonical_kmer, reverse_complement, SuperKmerInfos,
};

/// Retrieve the count of `superkmer` in `sk_count`.
/// the superkmer and its minimizer should be canonical.
/// Returns 0 if the superkmer is not found.
fn get_count_superkmer(sk_count: &SKCount, superkmer: &SuperKmerInfos) -> Count {
    let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
    for (hash, count) in sk_count.get_mut_iter(&superkmer.minimizer) {
        if *hash == superkmer_hash {
            return *count;
        }
    }
    // superkmer not found, count is 1
    0
}

/// Add 1 to the count of `superkmer` in `sk_count`.
/// If `superkmer` is not in `sk_count`, add it.
/// Returns the new count.
fn update_count_superkmer(sk_count: &mut SKCount, superkmer: &SuperKmerInfos) -> Count {
    let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
    for (super_kmer_hash, super_kmer_count) in sk_count.get_mut_iter(&superkmer.minimizer) {
        if *super_kmer_hash == superkmer_hash {
            let new_count = super_kmer_count.saturating_add(1);
            *super_kmer_count = new_count;
            return new_count;
        }
    }
    // superkmer not found, insert it, count is 1
    sk_count.insert(superkmer.minimizer.clone(), (superkmer_hash, 1));
    1
}

/// Searches for `hyperkmer_left` in `hk_count[minimizer]`
/// minimizer is assumed to be in canonical form
fn get_hyperkmer_left_id(
    hk_count: &HKCount,
    hyperkmers: &[String],
    minimizer: &str,
    hyperkmer_left: &str,
) -> Option<usize> {
    let hyperkmer_left = get_canonical_kmer(hyperkmer_left).1;
    for (left_hk, _right_hk, _count) in hk_count.get_iter(minimizer) {
        if hyperkmer_left == hyperkmers[left_hk.0] {
            return Some(left_hk.0);
        }
    }
    None
}

/// Searches for `hyperkmer_right` in `hk_count[minimizer]`
/// minimizer is assumed to be in canonical form
fn get_hyperkmer_right_id(
    hk_count: &HKCount,
    hyperkmers: &[String],
    minimizer: &str,
    hyperkmer_right: &str,
) -> Option<usize> {
    let hyperkmer_right = get_canonical_kmer(hyperkmer_right).1;
    for (_left_hk, right_hk, _count) in hk_count.get_iter(minimizer) {
        if hyperkmer_right == hyperkmers[right_hk.0] {
            return Some(right_hk.0);
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
    id_left_hk: (usize, bool),
    id_right_hk: (usize, bool),
    count: Count,
) {
    hk_count.insert(String::from(minimizer), (id_left_hk, id_right_hk, count));
}

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

fn get_left_and_rigth_from_hk(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> (String, String) {
    // Caution: the next and previous superkmer are given as the appear in the read.
    // To stay consistant accross all orientation of reading,
    // we conceptually switch the next and previous superkmers
    // if the current superkmer was not canonical in the read.

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
            assert!(current_sk.start_of_minimizer_as_read < previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read > next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk = current_sk.superkmer.len()
                - current_sk.start_of_minimizer_as_read
                + current_sk.start_of_superkmer_as_read
                - current_sk.minimizer.len();
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
    let end_of_left_hk = start_of_minimizer_in_sk + current_sk.minimizer.len() - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;
    let end_of_right_hk =
        start_of_minimizer_in_sk + distance_to_right + current_sk.minimizer.len() - 1;

    let left_hk = &current_sk.superkmer[start_of_left_hk..end_of_left_hk];
    let right_hk = &current_sk.superkmer[start_of_right_hk..end_of_right_hk];

    (String::from(left_hk), String::from(right_hk))
}

// TODO find a better name for the first stage function
fn first_stage(
    sequences: &Vec<&str>,
    k: usize,
    m: usize,
    threshold: Count,
) -> (SKCount, HKCount, Vec<String>) {
    let mut sk_count: SKCount = MashMap::new();
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

            // TODO stocker les counts pour ne pas les recalculer
            let current_count = update_count_superkmer(&mut sk_count, current_sk);

            // chain of comparisons ahead, but I don't need to be exhaustive and I find it good as it is
            // so I tell clippy to shup up
            #[allow(clippy::comparison_chain)]
            if current_count == threshold {
                let (left_hk, right_hk) =
                    get_left_and_rigth_from_hk(previous_sk, current_sk, next_sk);

                // left_hk and right_hk are the left and right hyperkmer
                // as we would see them if the minimizer was in canonical form in the read

                // OPTIMIZE maybe it is posssible to call get_hyperkmer_{left, right}_id and ignore get_count_superkmer
                // OPTIMIZE of even better: access the count of sk in streaming, so that no recomputation is needed
                let id_left_hk = if get_count_superkmer(&sk_count, previous_sk) >= threshold {
                    if previous_sk.was_read_canonical == current_sk.was_read_canonical {
                        get_hyperkmer_right_id(
                            &hk_count,
                            &hyperkmers,
                            &previous_sk.minimizer,
                            &left_hk,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        get_hyperkmer_left_id(
                            &hk_count,
                            &hyperkmers,
                            &previous_sk.minimizer,
                            &left_hk,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &left_hk)
                };
                let id_right_hk = if get_count_superkmer(&sk_count, next_sk) >= threshold {
                    if current_sk.was_read_canonical == next_sk.was_read_canonical {
                        get_hyperkmer_left_id(&hk_count, &hyperkmers, &next_sk.minimizer, &right_hk)
                            .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        get_hyperkmer_right_id(
                            &hk_count,
                            &hyperkmers,
                            &next_sk.minimizer,
                            &right_hk,
                        )
                        .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &right_hk)
                };

                // we have id left and rigth
                // let's get their orientation wrt to the orientation of the minimizer
                let left_hk = get_canonical_kmer(&left_hk);
                let right_hk = get_canonical_kmer(&right_hk);
                let left_change_orientation = left_hk.0 != current_sk.was_read_canonical;
                let right_change_orientation = right_hk.0 != current_sk.was_read_canonical;

                insert_new_entry_in_hyperkmer_count(
                    &mut hk_count,
                    &current_sk.minimizer,
                    (id_left_hk, left_change_orientation),
                    (id_right_hk, right_change_orientation),
                    current_count,
                );
            } else if current_count > threshold {
                let mut found = false;
                for (id_left_hk, id_right_hk, count_hk) in
                    hk_count.get_mut_iter(&current_sk.minimizer)
                {
                    let mut left_hk = hyperkmers[id_left_hk.0].clone();
                    if id_left_hk.1 {
                        left_hk = reverse_complement(&left_hk);
                    }

                    let mut right_hk = hyperkmers[id_right_hk.0].clone();
                    if id_right_hk.1 {
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
    sk_count: &mut SKCount,
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
            if get_count_superkmer(sk_count, superkmer) >= threshold {
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
            let left_hk_of_sk = &superkmer.superkmer[0..=(start_of_minimizer_in_sk + m - 1)];
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
    SKCount,
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

    let k = 100;
    let m = 31;
    let threshold = 2;

    let (_sk_count, _hk_count, _hyperkmers, _truncated_hk, _discarded_minimizers) =
        index_hyperkmers(k, m, threshold, &sequences);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_count_superkmer() {
        let sk0 = SuperKmerInfos {
            superkmer: "ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGG".into(),
            minimizer: "AATCGTCATCCGAAACGTCA".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 7,
            start_of_superkmer_as_read: 0,
        };
        let sk1 = SuperKmerInfos {
            superkmer: "GACGTTTCGGATGACGATTGTACGTGACGGTG".into(),
            minimizer: "ACAATCGTCATCCGAAACGT".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 9,
            start_of_superkmer_as_read: 8,
        };
        let sk2 = SuperKmerInfos {
            superkmer: "CGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATG".into(),
            minimizer: "ACCGTCACGTACAATCGTCA".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 19,
            start_of_superkmer_as_read: 10,
        };
        let sk3 = SuperKmerInfos {
            superkmer: "GACGATTGTACGTGACGGTGCGTCCGGATGAC".into(),
            minimizer: "ACGATTGTACGTGACGGTGC".into(),
            was_read_canonical: true,
            start_of_minimizer_as_read: 21,
            start_of_superkmer_as_read: 20,
        };

        let mut sk_count: SKCount = MashMap::new();

        for sk in vec![sk0, sk1, sk2, sk3] {
            assert_eq!(get_count_superkmer(&sk_count, &sk), 0);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 1);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 2);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 3);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 4);
            assert_eq!(get_count_superkmer(&sk_count, &sk), 4);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 5);
            assert_eq!(update_count_superkmer(&mut sk_count, &sk), 6);
            assert_eq!(get_count_superkmer(&sk_count, &sk), 6);
        }
    }

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
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(
                    minimizer.clone(),
                    ((random_left, true), (random_rigth, true), random_count),
                );
            }
        }

        // then link minimizer and hyperkmers
        for (i, minimizer) in minimizers.iter().enumerate() {
            let random_count: u16 = rng.gen_range(0..hyperkmers.len()) as u16;
            hk_count.insert(minimizer.clone(), ((i, true), (i + 1, true), random_count));
        }

        // add another random values
        for _ in 0..10 {
            for minimizer in &minimizers {
                let random_left = rng.gen_range(0..hyperkmers.len());
                let random_rigth = rng.gen_range(0..hyperkmers.len());
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(
                    minimizer.clone(),
                    ((random_left, true), (random_rigth, true), random_count),
                );
            }
        }

        // among all the random values, we are still able to get back our real data
        let hk_id = get_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizers[5], &hyperkmers[5]);
        assert_eq!(hk_id, Some(5));
        let hk_id = get_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizers[5], &hyperkmers[6]);
        assert_eq!(hk_id, Some(6));

        // query a non existant minimiser
        let minimizer: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(20)
            .map(char::from)
            .collect();
        let hk_id = get_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[5]);
        assert_eq!(hk_id, None);
        let hk_id = get_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[6]);
        assert_eq!(hk_id, None);

        // query an existing minimizer, but with a random hyperkmer
        let hyperkmer: String = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(100)
            .map(char::from)
            .collect();
        let hk_id = get_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmer);
        assert_eq!(hk_id, None);
        let hk_id = get_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmer);
        assert_eq!(hk_id, None);

        // query an existing minimizer, but with a random existing hyperkmer (TODO probability of failure is small but not null)
        let hk_id = get_hyperkmer_left_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[3]);
        assert_eq!(hk_id, None);
        let hk_id = get_hyperkmer_right_id(&hk_count, &hyperkmers, &minimizer, &hyperkmers[4]);
        assert_eq!(hk_id, None);
    }
}
