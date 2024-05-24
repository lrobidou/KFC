use log::warn;
use mashmap::MashMap;
use std::collections::HashMap;
use superkmers::{compute_superkmers, SuperKmerInfos};
use xxhash_rust::const_xxh3::xxh3_64;

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;

mod superkmers;

fn get_count_superkmer(
    sk_count: &MashMap<Minimizer, (HashSuperKmer, Count)>,
    superkmer: &SuperKmerInfos,
) -> Count {
    let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
    for (super_kmer_hash, super_kmer_count) in sk_count.get_mut_iter(&superkmer.minimizer) {
        if *super_kmer_hash == superkmer_hash {
            return *super_kmer_count;
        }
    }
    0
}

fn update_count_superkmer(
    sk_count: &mut MashMap<Minimizer, (HashSuperKmer, Count)>,
    superkmer: &SuperKmerInfos,
) -> Count {
    let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
    for (super_kmer_hash, super_kmer_count) in sk_count.get_mut_iter(&superkmer.minimizer) {
        if *super_kmer_hash == superkmer_hash {
            let new_count = *super_kmer_count + 1;
            *super_kmer_count = new_count;
            return new_count;
        }
    }
    // new_count is 1
    sk_count.insert(superkmer.minimizer.clone(), (superkmer_hash, 1));
    1
}

fn get_hyperkmer_left_id(
    hk_count: &MashMap<Minimizer, (usize, usize, Count)>,
    hyperkmers: &[String],
    minimizer: &str,
    hyperkmer_left: &str,
) -> Option<usize> {
    for (id_left_hk, _id_right_hk, _count) in hk_count.get_iter(&String::from(minimizer)) {
        if hyperkmer_left == hyperkmers[*id_left_hk] {
            return Some(*id_left_hk);
        }
    }
    None
}

fn get_hyperkmer_right_id(
    hk_count: &MashMap<Minimizer, (usize, usize, Count)>,
    hyperkmers: &[String],
    minimizer: &str,
    hyperkmer_right: &str,
) -> Option<usize> {
    for (_id_left_hk, id_right_hk, _count) in hk_count.get_iter(&String::from(minimizer)) {
        if hyperkmer_right == hyperkmers[*id_right_hk] {
            return Some(*id_right_hk);
        }
    }
    None
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

// TODO find a better name for the first stage function
fn first_stage(
    sequences: &Vec<&str>,
    k: usize,
    m: usize,
    threshold: Count,
) -> (
    MashMap<Minimizer, (HashSuperKmer, Count)>,
    MashMap<Minimizer, (usize, usize, Count)>,
    Vec<String>,
) {
    let mut sk_count: MashMap<Minimizer, (HashSuperKmer, Count)> = MashMap::new();
    let mut hk_count: MashMap<Minimizer, (usize, usize, Count)> = MashMap::new();
    let mut hyperkmers: Vec<String> = Vec::new();

    for sequence in sequences {
        let superkmers = compute_superkmers(sequence, k, m);
        for triplet in superkmers.windows(3) {
            // OPTIMIZE this implem is surely slow, we can do better by not using `String` everywhere, but &str
            let (previous_sk, current_sk, next_sk) = (&triplet[0], &triplet[1], &triplet[2]);

            // left and right parts of the current superkmer
            // truncated so that they do not include the current minimizer completely
            // (just like hyperkmers)
            let left_of_sk =
                &sequence[current_sk.start_of_superkmer..=current_sk.start_of_minimizer + m - 1];
            let right_of_sk = &sequence[current_sk.start_of_minimizer + 1
                ..current_sk.start_of_superkmer + current_sk.superkmer.len()];

            let current_count = update_count_superkmer(&mut sk_count, current_sk);

            // chain of comparisons ahead, but I don't need to be exhaustive and I find it good as it is
            // so I tell clippy to shup up
            #[allow(clippy::comparison_chain)]
            if current_count == threshold {
                let left_hk = {
                    let start_left_hk = previous_sk.start_of_minimizer + 1;
                    let end_left_hk = current_sk.start_of_minimizer + m - 1;
                    &sequence[start_left_hk..=end_left_hk]
                };
                let right_hk = {
                    let start_right_hk = current_sk.start_of_minimizer + 1;
                    let end_right_hk = next_sk.start_of_minimizer + m - 1;
                    &sequence[start_right_hk..=end_right_hk]
                };

                // OPTIMIZE maybe it is posssible to call get_hyperkmer_{left, right}_id and ignore get_count_superkmer
                let id_left_hk = if get_count_superkmer(&sk_count, previous_sk) >= threshold {
                    get_hyperkmer_right_id(&hk_count, &hyperkmers, &previous_sk.minimizer, left_hk)
                        .unwrap()
                } else {
                    let id_left_hk = hyperkmers.len();
                    hyperkmers.push(String::from(left_hk));
                    id_left_hk
                };
                let id_right_hk = if get_count_superkmer(&sk_count, next_sk) >= threshold {
                    get_hyperkmer_left_id(&hk_count, &hyperkmers, &current_sk.minimizer, right_hk)
                        .unwrap()
                } else {
                    let id_right_hk = hyperkmers.len();
                    hyperkmers.push(String::from(right_hk));
                    id_right_hk
                };
                hk_count.insert(
                    current_sk.minimizer.clone(),
                    (id_left_hk, id_right_hk, current_count),
                )
            } else if current_count > threshold {
                for (id_left_hk, id_right_hk, count_hk) in
                    hk_count.get_mut_iter(&current_sk.minimizer)
                {
                    let left_hk = &hyperkmers[*id_left_hk];
                    let right_hk = &hyperkmers[*id_right_hk];

                    if left_of_sk.ends_with(left_hk) && right_of_sk.starts_with(right_hk) {
                        *count_hk += 1;
                        break;
                    }
                }
            }
        }
    }
    (sk_count, hk_count, hyperkmers)
}

// TODO find a better name for the second stage function
fn second_stage(
    sk_count: &mut MashMap<Minimizer, (HashSuperKmer, Count)>,
    hk_count: &mut MashMap<Minimizer, (usize, usize, Count)>,
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
        let superkmers = compute_superkmers(sequence, k, m);
        for superkmer in &superkmers {
            if get_count_superkmer(sk_count, superkmer) >= threshold {
                continue;
            }
            let minimizer = &superkmer.minimizer;

            if !hk_count.contains_key(minimizer) {
                // increase count, set to 1 if it was 0
                let count = discarded_minimizers
                    .entry(minimizer.clone())
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
                if *count == threshold {
                    warn!(
                        "minimizer {} of superkmer {} is found {} times but its hyperkmer is not",
                        minimizer, superkmer.superkmer, count
                    );
                }
            }

            // left hyper-k-mer of sk
            // TODO I am taking more than the hyperkmer, but it should not be a problem since the match will stop at the hyperkmer anyway
            let start_of_minimizer_in_sk =
                superkmer.start_of_minimizer - superkmer.start_of_superkmer;
            let left_hk_of_sk = &superkmer.superkmer[0..=(start_of_minimizer_in_sk + m - 1)];
            let right_hk_of_sk =
                &superkmer.superkmer[(start_of_minimizer_in_sk + 1)..superkmer.superkmer.len()];
            let mut match_left = k - 1;
            let mut match_right = k - 1;
            let mut id_match_left = None;
            let mut id_match_right = None;

            for (id_left_hk, id_right_hk, _count) in hk_count.get_iter(minimizer) {
                let len_current_match_left =
                    common_suffix_length(left_hk_of_sk, &hyperkmers[*id_left_hk]);
                if len_current_match_left > match_left {
                    match_left = len_current_match_left;
                    id_match_left = Some(id_left_hk)
                }

                let len_current_match_right =
                    common_prefix_length(right_hk_of_sk, &hyperkmers[*id_right_hk]);
                if len_current_match_right > match_right {
                    match_right = len_current_match_right;
                    id_match_right = Some(id_right_hk)
                }
            }

            if let Some(id_match_left) = id_match_left {
                truncated_hk.insert(minimizer.clone(), (*id_match_left, match_left));
            }

            if let Some(id_match_right) = id_match_right {
                truncated_hk.insert(minimizer.clone(), (*id_match_right, match_left));
            }
        }
    }
    (truncated_hk, discarded_minimizers)
}

fn main() {
    let sequences =
        vec!["ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATGACGACGTACGTGACGGGGTCGGATGACG"];
    let k = 31;
    let m = 20;
    let threshold = 2;

    let (mut sk_count, mut hk_count, mut hyperkmers) = first_stage(&sequences, k, m, threshold);

    let (truncated_hk, discarded_minimizers) = second_stage(
        &mut sk_count,
        &mut hk_count,
        &hyperkmers,
        &sequences,
        k,
        m,
        threshold,
    );
    for sequence in sequences {
        let superkmers = superkmers::compute_superkmers(sequence, k, m);
        for superkmerinfo in superkmers {
            println!("{:?}", superkmerinfo);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_count_superkmer() {
        let sk0 = SuperKmerInfos {
            superkmer: "ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGG".into(),
            minimizer: "AATCGTCATCCGAAACGTCA".into(),
            is_forward: false,
            start_of_minimizer: 7,
            start_of_superkmer: 0,
        };
        let sk1 = SuperKmerInfos {
            superkmer: "GACGTTTCGGATGACGATTGTACGTGACGGTG".into(),
            minimizer: "ACAATCGTCATCCGAAACGT".into(),
            is_forward: false,
            start_of_minimizer: 9,
            start_of_superkmer: 8,
        };
        let sk2 = SuperKmerInfos {
            superkmer: "CGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATG".into(),
            minimizer: "ACCGTCACGTACAATCGTCA".into(),
            is_forward: false,
            start_of_minimizer: 19,
            start_of_superkmer: 10,
        };
        let sk3 = SuperKmerInfos {
            superkmer: "GACGATTGTACGTGACGGTGCGTCCGGATGAC".into(),
            minimizer: "ACGATTGTACGTGACGGTGC".into(),
            is_forward: true,
            start_of_minimizer: 21,
            start_of_superkmer: 20,
        };

        let mut sk_count: MashMap<Minimizer, (HashSuperKmer, Count)> = MashMap::new();

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
}
