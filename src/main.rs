use compute_left_and_right::{get_left_and_rigth_extended_hk, get_left_and_rigth_of_sk};
use fastxgz::fasta_reads;
use log::warn;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use std::{collections::HashMap, path::Path};

type Minimizer = String; // TODO change to integer when done
type HashSuperKmer = u64;
type Count = u16;

mod brrr_minimizers;
mod compute_left_and_right;
mod hyperkmers_counts;
mod search;
mod superkmer;
mod superkmers_computation;

use superkmer::{SubsequenceMetadata, Superkmer};

use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
use superkmers_computation::compute_superkmers_linear;

mod superkmers_count;

use superkmers_count::SuperKmerCounts;

fn extract_hyperkmer<'a>(
    hyperkmers: &'a [String],
    metadata: &HKMetadata,
) -> SubsequenceMetadata<'a> {
    let hyperkmer: &str = &hyperkmers[metadata.index];

    SubsequenceMetadata::new(
        hyperkmer,
        metadata.start,
        metadata.end,
        !metadata.change_orientation,
    )
}

/// Adds `new_hyperkmer` in `hyperkmers` and return its index
/// `new_hyperkmer` does not have to be in canonical form
fn add_new_hyperkmer(hyperkmers: &mut Vec<String>, new_hyperkmer: &SubsequenceMetadata) -> usize {
    let id_hyperkmer = hyperkmers.len();
    let new_hyperkmer_seq = new_hyperkmer.to_canonical_string();
    hyperkmers.push(new_hyperkmer_seq);
    id_hyperkmer
}

fn search_exact_ext_hyperkmer_match(
    hyperkmers: &[String],
    left_hk: &SubsequenceMetadata,
    right_hk: &SubsequenceMetadata,
    left_metadata: &HKMetadata,
    right_metadata: &HKMetadata,
) -> bool {
    // get sequences as they would appear if the current superkmer was canonical
    let candidate_left_hk = extract_hyperkmer(hyperkmers, left_metadata);
    let candidate_right_hk = extract_hyperkmer(hyperkmers, right_metadata);

    let match_left = candidate_left_hk.equal_str(&left_hk.to_string());
    let match_right = candidate_right_hk.equal_str(&right_hk.to_string());

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
    let mut hk_count: HKCount = HKCount::new();
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

            let (previous_sk, next_sk) = if current_sk.is_canonical_in_the_read() {
                (previous_sk, next_sk)
            } else {
                (next_sk, previous_sk)
            };
            // now, the beginning of current_sk.superkmer is close to the left neighbour

            // compute now if the previous and/or next superkmer is solid
            // (we do it here and now because the incoming instruction `increase_count_superkmer(current_sk)`
            // can increase their count as well if they are the same superkmer)
            let previous_sk_is_solid = sk_count.get_count_superkmer(previous_sk) >= threshold;
            let next_sk_is_solid = sk_count.get_count_superkmer(next_sk) >= threshold;

            // TODO est-il possible de stocker les counts pour ne pas les recalculer ?
            let current_count = sk_count.increase_count_superkmer(current_sk);

            // chain of comparisons ahead, but I don't need to be exhaustive and I find it good as it is
            // so I tell clippy to shup up
            #[allow(clippy::comparison_chain)]
            if current_count == threshold {
                let (left_extended_hk, right_extended_hk) =
                    get_left_and_rigth_extended_hk(previous_sk, current_sk, next_sk, k);

                // let canonical_extended_left_hk = left_extended_hk.0.to_canonical();
                // let canonical_extended_right_hk = right_extended_hk.0.to_canonical();

                // left_hk and right_hk are the left and right hyperkmer
                // as we would see them if the minimizer was in canonical form in the read

                // OPTIMIZE maybe it is posssible to call get_hyperkmer_{left, right}_id and ignore get_count_superkmer
                // OPTIMIZE of even better: access the count of sk in streaming, so that no recomputation is needed

                // OPTIMIZE the checks `previous_sk.was_read_canonical == current_sk.was_read_canonical` should be "avoidable" now
                let id_left_hk = if previous_sk_is_solid {
                    if previous_sk.is_canonical_in_the_read()
                        == current_sk.is_canonical_in_the_read()
                    {
                        hk_count
                            .get_extended_hyperkmer_right_id(
                                &hyperkmers,
                                &previous_sk.get_minimizer(),
                                &left_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        hk_count
                            .get_extended_hyperkmer_left_id(
                                &hyperkmers,
                                &previous_sk.get_minimizer(),
                                &left_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &left_extended_hk.0)
                };
                let id_right_hk = if next_sk_is_solid {
                    if current_sk.is_canonical_in_the_read() == next_sk.is_canonical_in_the_read() {
                        hk_count
                            .get_extended_hyperkmer_left_id(
                                &hyperkmers,
                                &next_sk.get_minimizer(),
                                &right_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        hk_count
                            .get_extended_hyperkmer_right_id(
                                &hyperkmers,
                                &next_sk.get_minimizer(),
                                &right_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    add_new_hyperkmer(&mut hyperkmers, &right_extended_hk.0)
                };

                // we have two ids (left and rigth) of extended hyperkmers containing our left and right hyperkemr
                // let's get their orientation wrt to the orientation of the minimizer

                let left_change_orientation = !left_extended_hk.0.is_canonical();
                let right_change_orientation = !right_extended_hk.0.is_canonical();

                let candidate_left_ext_hk =
                    SubsequenceMetadata::whole_string(&hyperkmers[id_left_hk])
                        .change_orientation_if(left_change_orientation);
                let candidate_right_ext_hk =
                    SubsequenceMetadata::whole_string(&hyperkmers[id_right_hk])
                        .change_orientation_if(right_change_orientation);

                assert!(candidate_left_ext_hk.equal_str(&left_extended_hk.0.to_string()));
                assert!(candidate_right_ext_hk.equal_str(&right_extended_hk.0.to_string()));
                assert!(hyperkmers[id_left_hk] == left_extended_hk.0.to_canonical_string());
                assert!(hyperkmers[id_right_hk] == right_extended_hk.0.to_canonical_string());

                hk_count.insert_new_entry_in_hyperkmer_count(
                    &current_sk.get_minimizer(),
                    &HKMetadata {
                        index: id_left_hk,
                        start: left_extended_hk.1,
                        end: left_extended_hk.2,
                        change_orientation: left_change_orientation,
                    },
                    &HKMetadata {
                        index: id_right_hk,
                        start: right_extended_hk.1,
                        end: right_extended_hk.2,
                        change_orientation: right_change_orientation,
                    },
                    current_count,
                );
            } else if current_count > threshold {
                let (left_sk, right_sk) = get_left_and_rigth_of_sk(current_sk);
                let found = hk_count.increase_count_if_exact_match(
                    current_sk,
                    &hyperkmers,
                    &left_sk,
                    &right_sk,
                );
                if !found {
                    // if no exact match, then we must at least have an approximate match
                    let new_left_and_right_metadata =
                        hk_count.search_for_inclusion(&hyperkmers, current_sk, &left_sk, &right_sk);
                    // TODO error message here
                    let (metadata_to_insert_left, metadata_to_insert_right) =
                        new_left_and_right_metadata.unwrap();
                    hk_count.insert_new_entry_in_hyperkmer_count(
                        &current_sk.get_minimizer(),
                        &metadata_to_insert_left,
                        &metadata_to_insert_right,
                        1,
                    );
                    // ensure the hyperkmer was correctly inserted
                    assert!(search_exact_hyperkmer_match(
                        &hyperkmers,
                        &left_sk,
                        &right_sk,
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
) -> HashMap<Minimizer, Count> {
    let mut discarded_minimizers: HashMap<Minimizer, Count> = HashMap::new();

    for sequence in sequences {
        let superkmers = compute_superkmers_linear(sequence, k, m);
        for superkmer in &superkmers {
            if sk_count.get_count_superkmer(superkmer) >= threshold {
                continue;
            }

            let minimizer = &superkmer.get_minimizer();

            if !hk_count.contains_minimizer(minimizer) {
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
                        minimizer,
                        superkmer.superkmer.to_string(),
                        count
                    );
                }
            }

            let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);
            let match_metadata = hk_count
                .search_for_maximal_inclusion(hyperkmers, k, minimizer, &left_sk, &right_sk);

            if let Some(metadata) = match_metadata {
                hk_count.insert_new_entry_in_hyperkmer_count(
                    minimizer,
                    &metadata.0,
                    &metadata.1,
                    1,
                );
            }
        }
    }
    discarded_minimizers
}

fn index_hyperkmers(
    k: usize,
    m: usize,
    threshold: Count,
    sequences: &Vec<&str>,
) -> (SuperKmerCounts, HKCount, Vec<String>, HashMap<String, u16>) {
    let start_fisrt_step = Instant::now();
    let (mut sk_count, mut hk_count, hyperkmers) = first_stage(sequences, k, m, threshold);
    println!(
        "time first stage: {} milliseconds",
        start_fisrt_step.elapsed().as_millis()
    );
    let start_second_stage = Instant::now();
    let discarded_minimizers = second_stage(
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
    (sk_count, hk_count, hyperkmers, discarded_minimizers)
}

fn stats_hk(hyperkmers: &[String], k: usize) -> (usize, usize) {
    for ext_hk in hyperkmers {
        assert_eq!(ext_hk.len(), k - 1);
    }
    (hyperkmers.len(), hyperkmers.len() * (k - 1))
}

fn compare_to_kmc(
    hk_count: &HKCount,
    hyperkmers: &[String],
    kmc_file: &Path,
    k: usize,
    m: usize,
    threshold: Count,
) {
    println!("mass comparison using kmc");
    let file = File::open(kmc_file).unwrap();
    let reader = BufReader::new(file);

    let (mut ok, mut ko) = (0, 0);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split = line.split('\t');
        let kmer = &split.next().unwrap();
        let count = &split.next().unwrap();
        let ground_truth_count: usize = count.parse().unwrap();

        if ground_truth_count >= threshold as usize {
            let kfc_res = search::search_kmer(hk_count, hyperkmers, kmer, k, m);
            if kfc_res as usize != ground_truth_count {
                println!("error on {kmer}, {kfc_res} != {ground_truth_count}");
                ko += 1;
            } else {
                ok += 1;
            }
        }
    }
    println!("ok = {ok}\nko={ko}");
}

fn dump_hk(lines: &[String]) {
    use std::io::Write;
    let mut file = File::create("hk.txt").unwrap();

    // Iterate over the vector and write each line to the file
    for line in lines {
        writeln!(file, "{line}").unwrap();
    }
}

fn main() {
    simple_logger::SimpleLogger::new().env().init().unwrap();
    let sequences: Vec<String> = fasta_reads("data/U00096.3.fasta")
        .unwrap()
        .map(|rcstring| rcstring.to_string())
        .collect();
    let sequences: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();

    let k = 100;
    let m = 20;
    let threshold = 1;

    let (_sk_count, hk_count, hyperkmers, _discarded_minimizers) =
        index_hyperkmers(k, m, threshold, &sequences);

    // let kmer_test = "CGCGAGGAGCTGGCCGAGGTGGATGTGGACTGGCTGATCGCCGAGCGCCCCGGCAAGGTAAGAACCTTGAAACAGCATCCACGCAAGAACAAAACGGCCA";

    // let count = search::search_kmer(&hk_count, &hyperkmers, kmer_test, k, m);
    // println!(
    //     "test with a kmer with a count 11 in KMC\nour count is {:?}",
    //     count
    // );
    dump_hk(&hyperkmers);
    let kmc_file_path = Path::new("100mers.txt");
    compare_to_kmc(&hk_count, &hyperkmers, kmc_file_path, k, m, threshold);
}

// TODO tests
// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_suffix_preffix() {
//         let x = "aazerty______";
//         let y = "aazerty-___---------";
//         assert_eq!(common_prefix_length(x, y), 7);
//         assert_eq!(common_prefix_length(y, x), 7);

//         assert_eq!(common_prefix_length("x", "y"), 0);
//         assert_eq!(common_prefix_length("", ""), 0);
//         assert_eq!(common_prefix_length("xyyyyyyyyyy", "xyy"), 3);

//         // suffix
//         let x = "aazerty_____-erc---------";
//         let y = "aazerty-___erc---------";
//         assert_eq!(common_suffix_length(x, y), 12);
//         assert_eq!(common_suffix_length(y, x), 12);

//         assert_eq!(common_suffix_length("x", "y"), 0);
//         assert_eq!(common_suffix_length("", ""), 0);
//         assert_eq!(common_suffix_length("xyyyyyyyyyy", "xyy"), 2);
//     }

//     #[test]
//     fn test_get_left_and_rigth_extended_hk() {
//         // in the sequence CAGATGGTTCAACCCTTAAGTTAGCGCTTATGGGATCAC
//         //     previous sk :              CCTTAAGTTAGCGCTTATGGGATCAC (GTTAGCGCTT)
//         //     current sk :            ACCCTTAAGTTAGCGCTTATG (CCCTTAAGTT)
//         //        next sk : CAGATGGTTCAACCCTTAAGTTAGCGCTTA (AACCCTTAAG)

//         //    SuperKmerInfos { superkmer: "GTGATCCCATAAGCGCTAACTTAAGG", minimizer: "AAGCGCTAAC", was_read_canonical: false, start_of_minimizer_as_read: 16715, start_of_superkmer_as_read: 16709 }
//         //    SuperKmerInfos { superkmer: "CATAAGCGCTAACTTAAGGGT", minimizer: "AACTTAAGGG", was_read_canonical: false, start_of_minimizer_as_read: 16708, start_of_superkmer_as_read: 16707 }
//         //    SuperKmerInfos { superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA", minimizer: "AACCCTTAAG", was_read_canonical: true, start_of_minimizer_as_read: 16706, start_of_superkmer_as_read: 16696 }
//         //    ("ACTTAAGGGT", "TAAGCGCTAACTTAAGGGT")
//         //    left extended hk: ("AGCGCTAACTTAAGG", 15)
//         //    right extended hk: ("TAAGCGCTAACTTAAGGGT", 10)
//         //    TAAGCGCTAACTTAAGGGTTGAACCATCTG
//         let previous_sk = SuperKmerInfos {
//             superkmer: "GTGATCCCATAAGCGCTAACTTAAGG".into(),
//             minimizer: "AAGCGCTAAC".into(),
//             was_read_canonical: false,
//             start_of_minimizer_as_read: 16715,
//             start_of_superkmer_as_read: 16709,
//         };
//         let current_sk = SuperKmerInfos {
//             superkmer: "CATAAGCGCTAACTTAAGGGT".into(),
//             minimizer: "AACTTAAGGG".into(),
//             was_read_canonical: false,
//             start_of_minimizer_as_read: 16708,
//             start_of_superkmer_as_read: 16707,
//         };
//         let next_sk = SuperKmerInfos {
//             superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA".into(),
//             minimizer: "AACCCTTAAG".into(),
//             was_read_canonical: true,
//             start_of_minimizer_as_read: 16706,
//             start_of_superkmer_as_read: 16696,
//         };
//         let (l, r) = get_left_and_rigth_extended_hk(&previous_sk, &current_sk, &next_sk);
//         //   left: "TAAGCGCTAACTTAAGGGT"
//         //  right: "ACCCTTAAGTTAGCGCTTA"
//         assert_eq!(r.0, "TAAGCGCTAACTTAAGGGT");
//         assert_eq!(l.0, "CATAAGCGCTAACTTAAGG");
//         //    sk right: ACCCTTAAGTTAGCGCTTA / TAAGCGCTAACTTAAGGGT
//         //    sk left: CCTTAAGTTAGCGCTTATG / CATAAGCGCTAACTTAAGG
//     }

//     #[test]
//     fn test_get_hyperkmer_left_and_rigth_id() {
//         use rand::{distributions::Alphanumeric, Rng}; // 0.8
//         let mut hk_count: HKCount = HKCount::new();
//         let mut hyperkmers: Vec<String> = Vec::new();

//         let nb_insertions = 10000;

//         // random insertions in hyperkmers
//         for _ in 0..nb_insertions {
//             let s: String = rand::thread_rng()
//                 .sample_iter(&Alphanumeric)
//                 .take(100)
//                 .map(char::from)
//                 .collect();
//             hyperkmers.push(get_canonical_kmer(&s).1);
//         }

//         // random minimizers
//         let mut minimizers = Vec::new();
//         for _ in 0..(hyperkmers.len() - 1) {
//             let minimizer: String = rand::thread_rng()
//                 .sample_iter(&Alphanumeric)
//                 .take(20)
//                 .map(char::from)
//                 .collect();
//             minimizers.push(get_canonical_kmer(&minimizer).1);
//         }

//         // link minimizer and hyperkmers
//         // start by inserting random values
//         let mut rng = rand::thread_rng();
//         for _ in 0..3 {
//             for minimizer in &minimizers {
//                 let random_left = rng.gen_range(0..hyperkmers.len());
//                 let random_left_start = rng.gen_range(0..hyperkmers[random_left].len() / 2);
//                 let random_left_end =
//                     rng.gen_range(hyperkmers[random_left].len() / 2..hyperkmers[random_left].len());

//                 let random_right = rng.gen_range(0..hyperkmers.len());
//                 let random_rigth_start = rng.gen_range(0..hyperkmers[random_right].len() / 2);
//                 let random_rigth_end = rng
//                     .gen_range(hyperkmers[random_right].len() / 2..hyperkmers[random_right].len());

//                 let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

//                 hk_count.insert_new_entry_in_hyperkmer_count(
//                     minimizer,
//                     (random_left, random_left_start, random_left_end, true),
//                     (random_right, random_rigth_start, random_rigth_end, true),
//                     random_count,
//                 );
//             }
//         }

//         // then link minimizer and hyperkmers
//         for (i, minimizer) in minimizers.iter().enumerate() {
//             let random_count: u16 = rng.gen_range(0..hyperkmers.len()) as u16;
//             let overlap = minimizer.len();
//             hk_count.insert_new_entry_in_hyperkmer_count(
//                 minimizer,
//                 (i, 0, overlap, true),
//                 (i + 1, 0, overlap, true),
//                 random_count,
//             );
//         }

//         // add another random values
//         for _ in 0..10 {
//             for minimizer in &minimizers {
//                 let random_left = rng.gen_range(0..hyperkmers.len());
//                 let random_left_start = rng.gen_range(0..hyperkmers[random_left].len() / 2);
//                 let random_left_end =
//                     rng.gen_range(hyperkmers[random_left].len() / 2..hyperkmers[random_left].len());

//                 let random_right = rng.gen_range(0..hyperkmers.len());
//                 let random_rigth_start = rng.gen_range(0..hyperkmers[random_right].len() / 2);
//                 let random_rigth_end = rng
//                     .gen_range(hyperkmers[random_right].len() / 2..hyperkmers[random_right].len());

//                 let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

//                 hk_count.insert_new_entry_in_hyperkmer_count(
//                     minimizer,
//                     (random_left, random_left_start, random_left_end, true),
//                     (random_right, random_rigth_start, random_rigth_end, true),
//                     random_count,
//                 );
//             }
//         }

//         // among all the random values, we are still able to get back our real data
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizers[5], &hyperkmers[5]);
//         assert_eq!(hk_id, Some(5));
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizers[5], &hyperkmers[6]);
//         assert_eq!(hk_id, Some(6));

//         // query a non existant minimiser
//         let minimizer: String = rand::thread_rng()
//             .sample_iter(&Alphanumeric)
//             .take(20)
//             .map(char::from)
//             .collect();
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmers[5]);
//         assert_eq!(hk_id, None);
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmers[6]);
//         assert_eq!(hk_id, None);

//         // query an existing minimizer, but with a random hyperkmer
//         let hyperkmer: String = rand::thread_rng()
//             .sample_iter(&Alphanumeric)
//             .take(100)
//             .map(char::from)
//             .collect();
//         let hk_id = hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmer);
//         assert_eq!(hk_id, None);
//         let hk_id = hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmer);
//         assert_eq!(hk_id, None);

//         // query an existing minimizer, but with a random existing hyperkmer (TODO probability of failure is small but not null)
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmers[3]);
//         assert_eq!(hk_id, None);
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmers[4]);
//         assert_eq!(hk_id, None);
//     }

//     // #[test]
//     // fn test_suf_of_x_is_pref_of_y() {
//     //     assert_eq!(suf_of_x_is_pref_of_y("abcdefghij", "fghijyhy"), 5);
//     // }
// }
