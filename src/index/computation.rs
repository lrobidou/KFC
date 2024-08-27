use crate::{prt, superkmer::Superkmer};
use std::{
    collections::HashMap,
    path::Path,
    sync::{Arc, Mutex, RwLock},
};

use bitvec::array::IntoIter;
use fastxgz::fasta_reads;
use itertools::Itertools;
use log::warn;

use crate::{
    compute_left_and_right::{get_left_and_rigth_extended_hk, get_left_and_rigth_of_sk},
    index::components::{
        add_new_large_hyperkmer, get_subsequence_from_metadata, search_exact_hyperkmer_match,
        HKMetadata,
    },
    superkmers_computation::compute_superkmers_linear_streaming,
    Count, Minimizer,
};

use super::{
    components::{HKCount, ParallelExtendedHyperkmers, SuperKmerCounts},
    parallel::{
        // mac::{read_lock, write_lock},
        Parallel,
    },
};

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as likely;
#[cfg(feature = "nightly")]
use core::intrinsics::likely;

fn is_sk_solid(sk_count: &Parallel<SuperKmerCounts>, sk: &Superkmer, threshold: Count) -> bool {
    let sk_count = sk_count.get_from_minimizer(sk.get_minimizer());
    let sk_count = sk_count.read().unwrap();
    sk_count.get_count_superkmer(sk) >= threshold
}

pub fn first_stage<P: AsRef<Path>>(
    path: P,
    k: usize,
    m: usize,
    threshold: Count,
) -> (
    Parallel<SuperKmerCounts>,
    Parallel<HKCount>,
    Arc<RwLock<ParallelExtendedHyperkmers>>,
    Arc<RwLock<Vec<(usize, Vec<u8>)>>>,
) {
    let sk_count = Parallel::<SuperKmerCounts>::new(SuperKmerCounts::new);
    let hk_count = Parallel::<HKCount>::new(HKCount::new);
    let hyperkmers = Arc::new(RwLock::new(ParallelExtendedHyperkmers::new(k, 1000)));
    let large_hyperkmers: Arc<RwLock<Vec<(usize, Vec<u8>)>>> = Arc::new(RwLock::new(Vec::new()));

    // TODO change and iterate over u8, and also not a vec
    let sequences: Vec<String> = fasta_reads(path)
        .unwrap()
        .map(|rcstring| rcstring.to_string())
        .collect();

    rayon::scope(|s| {
        let chunks = sequences.into_iter().chunks(100);
        for chunk in chunks.into_iter() {
            let sk_count = sk_count.clone();
            let hk_count = hk_count.clone();
            let hyperkmers = hyperkmers.clone();
            let large_hyperkmers = large_hyperkmers.clone();
            let lines: Vec<String> = chunk.into_iter().collect_vec();
            let lines = Arc::new(Mutex::new(lines.clone()));
            s.spawn(move |_| {
                first_stage_for_a_chunck(
                    lines,
                    k,
                    m,
                    threshold,
                    &sk_count,
                    &hk_count,
                    &hyperkmers,
                    &large_hyperkmers,
                )
            });
        }
    });
    (sk_count, hk_count, hyperkmers, large_hyperkmers)
}

pub fn second_stage<P: AsRef<Path>>(
    sk_count: &mut Parallel<SuperKmerCounts>,
    hk_count: &mut Parallel<HKCount>,
    hyperkmers: Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: Arc<RwLock<Vec<(usize, Vec<u8>)>>>,
    path: P, // OPTIMIZE prendre un iterateur sur des &[u8] ?
    k: usize,
    m: usize,
    threshold: Count,
) -> Parallel<HashMap<u64, u16>> {
    let discarded_minimizers: Parallel<HashMap<u64, u16>> =
        Parallel::<HashMap<Minimizer, Count>>::new(HashMap::new);

    let sequences: Vec<String> = fasta_reads(path)
        .unwrap()
        .map(|rcstring| rcstring.to_string())
        .collect();

    rayon::scope(|s| {
        let chunks = sequences.into_iter().chunks(100);
        for chunk in chunks.into_iter() {
            let sk_count = sk_count.clone();
            let hk_count = hk_count.clone();
            let hyperkmers = hyperkmers.clone();
            let large_hyperkmers = large_hyperkmers.clone();
            let discarded_minimizers = discarded_minimizers.clone();
            let lines: Vec<String> = chunk.into_iter().collect_vec();
            let lines = Arc::new(Mutex::new(lines.clone()));
            s.spawn(move |_| {
                second_stage_for_a_chunk(
                    &sk_count,
                    &hk_count,
                    &hyperkmers,
                    &large_hyperkmers,
                    lines,
                    k,
                    m,
                    threshold,
                    &discarded_minimizers,
                )
            });
        }
    });
    discarded_minimizers
}

// TODO "style" find a better name for the first stage function
pub fn first_stage_for_a_chunck(
    sequences: Arc<Mutex<Vec<String>>>,
    k: usize,
    m: usize,
    threshold: Count,
    sk_count: &Parallel<SuperKmerCounts>,
    hk_count: &Parallel<HKCount>,
    hyperkmers: &Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: &Arc<RwLock<Vec<(usize, Vec<u8>)>>>,
) {
    let sequences = sequences.lock().unwrap();
    for sequence in sequences.iter() {
        let sequence = &sequence.as_bytes();
        let superkmers = match compute_superkmers_linear_streaming(sequence, k, m) {
            Some(superkmers_iter) => superkmers_iter,
            None => continue,
        };

        for (previous_sk, current_sk, next_sk) in superkmers.into_iter().tuple_windows() {
            let (previous_sk, next_sk) = if current_sk.is_canonical_in_the_read() {
                (previous_sk, next_sk)
            } else {
                (next_sk, previous_sk)
            };
            // now, the beginning of current_sk.superkmer is close to the left neighbour

            // compute now if the previous and/or next superkmer is solid
            // (we do it here and now because the incoming instruction `increase_count_superkmer(current_sk)`
            // can increase their count as well if they are the same superkmer)
            let previous_sk_is_solid = is_sk_solid(sk_count, &previous_sk, threshold);
            let next_sk_is_solid = is_sk_solid(sk_count, &next_sk, threshold);

            let current_sk_count = sk_count.get_from_minimizer(current_sk.get_minimizer());
            let mut current_sk_count = current_sk_count.write().unwrap();
            // OPTIMIZE est-il possible de stocker les counts pour ne pas les recalculer ?
            let current_count = current_sk_count.increase_count_superkmer(&current_sk);
            // chain of comparisons ahead, but I don't need to be exhaustive and I find it good as it is
            // so I tell clippy to shup up
            #[allow(clippy::comparison_chain)]
            if current_count == threshold {
                let (left_extended_hk, right_extended_hk) =
                    get_left_and_rigth_extended_hk(&previous_sk, &current_sk, &next_sk, k);

                // left_hk and right_hk are the left and right hyperkmer
                // as we would see them if the minimizer was in canonical form in the read

                // OPTIMIZE maybe it is posssible to call get_hyperkmer_{left, right}_id and ignore get_count_superkmer
                // OPTIMIZE of even better: access the count of sk in streaming, so that no recomputation is needed
                let (id_left_hk, is_large_left) = if previous_sk_is_solid {
                    let previous_hk_count =
                        hk_count.get_from_minimizer(previous_sk.get_minimizer());
                    let previous_hk_count = previous_hk_count.read().unwrap();
                    if previous_sk.is_canonical_in_the_read()
                        == current_sk.is_canonical_in_the_read()
                    {
                        previous_hk_count
                            .get_extended_hyperkmer_right_id(
                                &hyperkmers.read().unwrap(),
                                &large_hyperkmers.read().unwrap(),
                                &previous_sk.get_minimizer(),
                                &left_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        previous_hk_count
                            .get_extended_hyperkmer_left_id(
                                &hyperkmers.read().unwrap(),
                                &large_hyperkmers.read().unwrap(),
                                &previous_sk.get_minimizer(),
                                &left_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    // previous sk is not solid => our hyperkmer is not already present
                    // let's add it
                    if likely(!left_extended_hk.3) {
                        (
                            hyperkmers
                                .read()
                                .unwrap()
                                .add_new_ext_hyperkmer(&left_extended_hk.0),
                            false,
                        )
                    } else {
                        (
                            add_new_large_hyperkmer(
                                &mut large_hyperkmers.write().unwrap(),
                                &left_extended_hk.0,
                            ),
                            true,
                        )
                    }
                };
                let (id_right_hk, is_large_right) = if next_sk_is_solid {
                    let next_hk_count = hk_count.get_from_minimizer(next_sk.get_minimizer());
                    let next_hk_count = next_hk_count.read().unwrap();
                    if current_sk.is_canonical_in_the_read() == next_sk.is_canonical_in_the_read() {
                        next_hk_count
                            .get_extended_hyperkmer_left_id(
                                &hyperkmers.read().unwrap(),
                                &large_hyperkmers.read().unwrap(),
                                &next_sk.get_minimizer(),
                                &right_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    } else {
                        next_hk_count
                            .get_extended_hyperkmer_right_id(
                                &hyperkmers.read().unwrap(),
                                &large_hyperkmers.read().unwrap(),
                                &next_sk.get_minimizer(),
                                &right_extended_hk.0,
                            )
                            .expect("Hash collision on superkmers. Please change your seed.")
                    }
                } else {
                    // previous sk is not solid => our hyperkmer is not already present
                    // let's add it
                    if likely(!right_extended_hk.3) {
                        (
                            hyperkmers
                                .read()
                                .unwrap()
                                .add_new_ext_hyperkmer(&right_extended_hk.0),
                            false,
                        )
                    } else {
                        (
                            add_new_large_hyperkmer(
                                &mut large_hyperkmers.write().unwrap(),
                                &right_extended_hk.0,
                            ),
                            true,
                        )
                    }
                };

                // we have two ids (left and rigth) of extended hyperkmers containing our left and right hyperkemr
                // let's get their orientation wrt to the orientation of the minimizer

                let left_change_orientation = !left_extended_hk.0.is_canonical();
                let right_change_orientation = !right_extended_hk.0.is_canonical();

                let left_hk_metadata = HKMetadata::new(
                    id_left_hk,
                    left_extended_hk.1,
                    left_extended_hk.2,
                    is_large_left,
                    left_change_orientation,
                );

                let right_hk_metadata = HKMetadata::new(
                    id_right_hk,
                    right_extended_hk.1,
                    right_extended_hk.2,
                    is_large_right,
                    right_change_orientation,
                );
                let hyperkmers = hyperkmers.read().unwrap();
                let large_hyperkmers = large_hyperkmers.read().unwrap();
                let candidate_left_ext_hk = &get_subsequence_from_metadata(
                    &hyperkmers,
                    &large_hyperkmers,
                    &left_hk_metadata,
                )
                .change_orientation_if(left_change_orientation);
                let candidate_right_ext_hk = &get_subsequence_from_metadata(
                    &hyperkmers,
                    &large_hyperkmers,
                    &right_hk_metadata,
                )
                .change_orientation_if(right_change_orientation);

                debug_assert!(left_extended_hk.0.equal_bitpacked(candidate_left_ext_hk));
                debug_assert!(right_extended_hk.0.equal_bitpacked(candidate_right_ext_hk));
                debug_assert!(left_extended_hk.0.to_canonical().equal_bitpacked(
                    &get_subsequence_from_metadata(
                        &hyperkmers,
                        &large_hyperkmers,
                        &left_hk_metadata,
                    )
                ));
                debug_assert!(right_extended_hk.0.to_canonical().equal_bitpacked(
                    &get_subsequence_from_metadata(
                        &hyperkmers,
                        &large_hyperkmers,
                        &right_hk_metadata,
                    )
                ));

                #[cfg(debug_assertions)]
                {
                    let left_ext_hk = get_subsequence_from_metadata(
                        &hyperkmers,
                        &large_hyperkmers,
                        &left_hk_metadata,
                    )
                    .change_orientation_if(left_hk_metadata.get_change_orientation());

                    let right_ext_hk = get_subsequence_from_metadata(
                        &hyperkmers,
                        &large_hyperkmers,
                        &right_hk_metadata,
                    )
                    .change_orientation_if(right_hk_metadata.get_change_orientation());

                    // extract candidate hyperkmers
                    let left_hyperkmer = &left_ext_hk
                        .subsequence(left_hk_metadata.get_start(), left_hk_metadata.get_end());
                    let right_hyperkmer = &right_ext_hk
                        .subsequence(right_hk_metadata.get_start(), right_hk_metadata.get_end());

                    let left_string = left_hyperkmer.to_string();
                    let right_string = right_hyperkmer.to_string();

                    debug_assert_eq!(
                        left_string[(left_string.len() - (m - 2))..left_string.len()],
                        right_string[0..(m - 2)]
                    );
                }
                let current_hk_count = hk_count.get_from_minimizer(current_sk.get_minimizer());
                let mut current_hk_count = current_hk_count.write().unwrap();
                current_hk_count.insert_new_entry_in_hyperkmer_count(
                    &current_sk.get_minimizer(),
                    &left_hk_metadata,
                    &right_hk_metadata,
                    current_count,
                );
            } else if current_count > threshold {
                let current_hk_count = hk_count.get_from_minimizer(current_sk.get_minimizer());
                let mut current_hk_count = current_hk_count.write().unwrap();
                // TODO fusionnner les deux passes
                let hyperkmers = hyperkmers.read().unwrap();
                let large_hyperkmers = large_hyperkmers.read().unwrap();
                let (left_sk, right_sk) = get_left_and_rigth_of_sk(&current_sk);
                let found = current_hk_count.increase_count_if_exact_match(
                    &current_sk.get_minimizer(),
                    &hyperkmers,
                    &large_hyperkmers,
                    &left_sk,
                    &right_sk,
                );
                if !found {
                    // if no exact match, then we must at least have an approximate match
                    let new_left_and_right_metadata = current_hk_count.search_for_inclusion(
                        &hyperkmers,
                        &large_hyperkmers,
                        &current_sk,
                        &left_sk,
                        &right_sk,
                    );

                    // If we are here, the superkmer is solid. Therefore, it must have been inserted.
                    let (metadata_to_insert_left, metadata_to_insert_right) =
                        new_left_and_right_metadata
                            .expect("Hash collision on superkmers. Please change your seed.");
                    current_hk_count.insert_new_entry_in_hyperkmer_count(
                        &current_sk.get_minimizer(),
                        &metadata_to_insert_left,
                        &metadata_to_insert_right,
                        1,
                    );
                    // ensure the hyperkmer was correctly inserted
                    debug_assert!(search_exact_hyperkmer_match(
                        &hyperkmers,
                        &large_hyperkmers,
                        &left_sk,
                        &right_sk,
                        &metadata_to_insert_left,
                        &metadata_to_insert_right
                    ));
                }
            }
        }
    }
    // (sk_count, hk_count, hyperkmers, large_hyperkmers) // TODO "style" renommer extended_hyperkmers
}

// TODO "style" find a better name for the second stage function
pub fn second_stage_for_a_chunk(
    sk_count: &Parallel<SuperKmerCounts>,
    hk_count: &Parallel<HKCount>,
    hyperkmers: &Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: &Arc<RwLock<Vec<(usize, Vec<u8>)>>>,
    sequences: Arc<Mutex<Vec<String>>>, // OPTIMIZE prendre un iterateur sur des &[u8] ?
    k: usize,
    m: usize,
    threshold: Count,
    discarded_minimizers: &Parallel<HashMap<Minimizer, Count>>,
) {
    let hyperkmers = hyperkmers.read().unwrap();
    let large_hyperkmers = large_hyperkmers.read().unwrap();
    let sequences = sequences.lock().unwrap();
    for sequence in sequences.iter() {
        let sequence = &sequence.as_bytes();
        let superkmers = match compute_superkmers_linear_streaming(sequence, k, m) {
            Some(superkmers_iter) => superkmers_iter,
            None => continue,
        };
        for superkmer in superkmers {
            let minimizer = superkmer.get_minimizer();
            let sk_count = sk_count.get_from_minimizer(minimizer);
            let sk_count = sk_count.read().unwrap();

            if sk_count.get_count_superkmer(&superkmer) >= threshold {
                continue;
            }

            let hk_count = hk_count.get_from_minimizer(minimizer);
            let mut hk_count = hk_count.write().unwrap();

            let discarded_minimizers = discarded_minimizers.get_from_minimizer(minimizer);
            let mut discarded_minimizers = discarded_minimizers.write().unwrap();
            if !hk_count.contains_minimizer(&minimizer) {
                // increase count, set to 1 if it was 0
                let count = discarded_minimizers
                    .entry(minimizer)
                    .and_modify(|counter| {
                        *counter = counter.saturating_add(1);
                    })
                    .or_insert(1);
                if *count == threshold {
                    // TODO debug: what if t == 1 ? Then we miss the first and last superkmer
                    warn!(
                        "minimizer {} of superkmer {} is found {} times but its hyperkmer is not",
                        minimizer,
                        superkmer.superkmer.to_string(),
                        count
                    );
                }
            }

            let (left_sk, right_sk) = get_left_and_rigth_of_sk(&superkmer);
            let match_metadata = hk_count.search_for_maximal_inclusion(
                &hyperkmers,
                &large_hyperkmers,
                k,
                m,
                &minimizer,
                &left_sk,
                &right_sk,
            );

            // TODO duplication possible
            if let Some(metadata) = match_metadata {
                hk_count.insert_new_entry_in_hyperkmer_count(
                    &minimizer,
                    &metadata.0,
                    &metadata.1,
                    1,
                );
            }
        }
    }
}
