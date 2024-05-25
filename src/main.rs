use log::warn;
use mashmap::MashMap;
use std::collections::HashMap;
use std::time::Instant;
use xxhash_rust::const_xxh3::xxh3_64;

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;

mod superkmers;

use fastxgz::fasta_reads;
use superkmers::{compute_superkmers_precompute_mmers, SuperKmerInfos};

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
        let start_superkmers = Instant::now();
        let superkmers = compute_superkmers_precompute_mmers(sequence, k, m);
        println!(
            "super kmers computed in {} seconds",
            start_superkmers.elapsed().as_secs()
        );

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
                    get_hyperkmer_left_id(&hk_count, &hyperkmers, &next_sk.minimizer, right_hk)
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
        let superkmers = compute_superkmers_precompute_mmers(sequence, k, m); // TODO discuss if we compute them twice or store them accross steps
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
    let sequences: Vec<String> = fasta_reads("data/U00096.3.fasta")
        .unwrap()
        .map(|rcstring| rcstring.to_string())
        .collect();
    let sequences: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();

    // println!("{:?}", sequences);

    let k = 40;
    let m = 31;
    let threshold = 2;

    let start_fisrt_step = Instant::now();
    let (mut sk_count, mut hk_count, hyperkmers) = first_stage(&sequences, k, m, threshold);
    println!(
        "time first step: {} seconds",
        start_fisrt_step.elapsed().as_secs()
    );
    let start_second_stage = Instant::now();
    let (truncated_hk, discarded_minimizers) = second_stage(
        &mut sk_count,
        &mut hk_count,
        &hyperkmers,
        &sequences,
        k,
        m,
        threshold,
    );
    println!(
        "time second stage: {} seconds",
        start_second_stage.elapsed().as_secs()
    );
    // for sequence in sequences {
    //     let superkmers = superkmers::compute_superkmers_linear(sequence, k, m);
    //     for superkmerinfo in superkmers {
    //         println!("{:?}", superkmerinfo);
    //     }
    // }
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
        let mut hk_count: MashMap<Minimizer, (usize, usize, Count)> = MashMap::new();
        let mut hyperkmers: Vec<String> = Vec::new();

        let nb_insertions = 10000;

        // random insertions in hyperkmers
        for _ in 0..nb_insertions {
            let s: String = rand::thread_rng()
                .sample_iter(&Alphanumeric)
                .take(100)
                .map(char::from)
                .collect();
            hyperkmers.push(s);
        }

        // random minimizers
        let mut minimizers = Vec::new();
        for _ in 0..(hyperkmers.len() - 1) {
            let minimizer: String = rand::thread_rng()
                .sample_iter(&Alphanumeric)
                .take(20)
                .map(char::from)
                .collect();
            minimizers.push(minimizer);
        }

        // link minimizer and hyperkmers
        // start by inserting random values
        let mut rng = rand::thread_rng();
        for _ in 0..3 {
            for minimizer in &minimizers {
                let random_left = rng.gen_range(0..hyperkmers.len());
                let random_rigth = rng.gen_range(0..hyperkmers.len());
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(minimizer.clone(), (random_left, random_rigth, random_count));
            }
        }

        // then link minimizer and hyperkmers
        for (i, minimizer) in minimizers.iter().enumerate() {
            let random_count: u16 = rng.gen_range(0..hyperkmers.len()) as u16;
            hk_count.insert(minimizer.clone(), (i, i + 1, random_count));
        }

        // add another random values
        for _ in 0..10 {
            for minimizer in &minimizers {
                let random_left = rng.gen_range(0..hyperkmers.len());
                let random_rigth = rng.gen_range(0..hyperkmers.len());
                let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

                hk_count.insert(minimizer.clone(), (random_left, random_rigth, random_count));
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
