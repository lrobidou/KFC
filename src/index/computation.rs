use crate::{
    buckets::{Buckets, LockPosition, NB_BUCKETS},
    compute_left_and_right::{get_left_and_rigth_extended_hk, get_left_and_rigth_of_sk},
    index::{
        components::{
            add_new_hyperkmer, search_exact_hyperkmer_match, HKCount, HKMetadata,
            ParallelExtendedHyperkmers, SuperKmerCounts,
        },
        LargeExtendedHyperkmer,
    },
    subsequence::{NoBitPacked, Subsequence},
    superkmer::Superkmer,
    superkmers_computation::compute_superkmers_linear_streaming,
    Count, Minimizer,
};

use std::{
    collections::HashMap,
    path::Path,
    sync::{Arc, Mutex, RwLock, RwLockReadGuard, RwLockWriteGuard},
};

use itertools::Itertools;
use log::warn;
use needletail::{parse_fastx_file, FastxReader};

const BATCH_SIZE: usize = 100;

fn replace_n(sequences: &mut Vec<Vec<u8>>) {
    for sequence in sequences {
        for char in sequence {
            if char == &b'N' {
                *char = b'A';
            }
        }
    }
}

/// First stage of the construction of the KFC index.
///
/// Splits each line of the input into chunks, then performs the first stage on each chunk in parallel.  
pub fn first_stage<P: AsRef<Path>>(
    path: P,
    k: usize,
    m: usize,
    threshold: Count,
) -> (
    Buckets<SuperKmerCounts>,
    Buckets<HKCount>,
    Arc<RwLock<ParallelExtendedHyperkmers>>,
    Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
) {
    let sk_count = Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new);
    let hk_count = Buckets::<HKCount>::new(HKCount::new);
    let hyperkmers = Arc::new(RwLock::new(ParallelExtendedHyperkmers::new(k, 1000)));
    let large_hyperkmers: Arc<RwLock<Vec<LargeExtendedHyperkmer>>> =
        Arc::new(RwLock::new(Vec::new()));

    let sequences = parse_fastx_file(path).unwrap();
    let sequences = LinesIter::new(sequences);

    rayon::scope(|s| {
        let chunks = sequences.into_iter().chunks(BATCH_SIZE);
        for chunk in chunks.into_iter() {
            let sk_count = sk_count.clone();
            let hk_count = hk_count.clone();
            let hyperkmers = hyperkmers.clone();
            let large_hyperkmers = large_hyperkmers.clone();
            let lines = chunk.into_iter().collect_vec(); //TODO copy
            let lines = Arc::new(Mutex::new(lines));
            s.spawn(move |_| {
                let mut sequences = lines.lock().unwrap();
                first_stage_for_a_chunck(
                    &mut sequences,
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

/// Second stage of the construction of the KFC index.
///
/// Splits each line of the input into chunks, then performs the second stage on each chunk in parallel.  
#[allow(clippy::too_many_arguments)]
pub fn second_stage<P: AsRef<Path>>(
    sk_count: &Buckets<SuperKmerCounts>,
    hk_count: &Buckets<HKCount>,
    hyperkmers: Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
    path: P,
    k: usize,
    m: usize,
    threshold: Count,
) -> Buckets<HashMap<u64, u16>> {
    let discarded_minimizers = Buckets::<HashMap<Minimizer, Count>>::new(HashMap::new);

    let sequences = parse_fastx_file(path).unwrap();
    let sequences = LinesIter::new(sequences);

    rayon::scope(|s| {
        let chunks = sequences.into_iter().chunks(BATCH_SIZE);
        for chunk in chunks.into_iter() {
            let sk_count = sk_count.clone();
            let hk_count = hk_count.clone();
            let hyperkmers = hyperkmers.clone();
            let large_hyperkmers = large_hyperkmers.clone();
            let discarded_minimizers = discarded_minimizers.clone();
            let lines = chunk.into_iter().collect_vec(); // TODO copy
            let lines = Arc::new(Mutex::new(lines));

            s.spawn(move |_| {
                let mut sequences = lines.lock().unwrap();
                second_stage_for_a_chunk(
                    &sk_count,
                    &hk_count,
                    &hyperkmers,
                    &large_hyperkmers,
                    &mut sequences,
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
struct LinesIter {
    data: Box<dyn FastxReader>,
}

impl LinesIter {
    fn new(data: Box<dyn FastxReader>) -> Self {
        Self { data }
    }
}

impl Iterator for LinesIter {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.data.next()?.unwrap();
        let sequence = record.seq();
        let sequence_vec = sequence.iter().copied().collect(); // TODO copy
        Some(sequence_vec)
    }
}

fn get_bucket_of_previous_hk(
    current_sk: &Superkmer,
    previous_sk: &Superkmer,
    hk_count_locks: &(
        RwLockWriteGuard<HKCount>,
        Option<RwLockReadGuard<HKCount>>,
        Option<RwLockReadGuard<HKCount>>,
        LockPosition,
        LockPosition,
    ),
    left_extended_hk: &Subsequence<NoBitPacked>,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
) -> (usize, usize, bool) {
    // used to access buckets and compute hyperkmer id
    // read hk_count table associated with the previous minimizer
    let large_hyperkmers = large_hyperkmers.read().unwrap();
    let previous_minimizer = previous_sk.get_minimizer();

    let previous_hk_count: &HKCount = match hk_count_locks.3 {
        LockPosition::ThisLock => hk_count_locks.1.as_ref().unwrap(),
        LockPosition::CurrentLock => &hk_count_locks.0,
        LockPosition::OtherLock => hk_count_locks.2.as_ref().unwrap(),
    };

    if previous_sk.is_canonical_in_the_read() == current_sk.is_canonical_in_the_read() {
        previous_hk_count
            .get_extended_hyperkmer_right_id(
                hyperkmers,
                &large_hyperkmers,
                &previous_minimizer,
                left_extended_hk,
            )
            .expect("Hash collision on superkmers. Please change your seed.")
    } else {
        previous_hk_count
            .get_extended_hyperkmer_left_id(
                hyperkmers,
                &large_hyperkmers,
                &previous_minimizer,
                left_extended_hk,
            )
            .expect("Hash collision on superkmers. Please change your seed.")
    }
}

fn get_bucket_of_next_hk(
    current_sk: &Superkmer,
    next_sk: &Superkmer,
    hk_count_locks: &(
        RwLockWriteGuard<HKCount>,
        Option<RwLockReadGuard<HKCount>>,
        Option<RwLockReadGuard<HKCount>>,
        LockPosition,
        LockPosition,
    ),
    right_extended_hk: &Subsequence<NoBitPacked>,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
) -> (usize, usize, bool) {
    let next_minimizer = next_sk.get_minimizer();
    let large_hyperkmers = large_hyperkmers.read().unwrap();
    let next_hk_count: &HKCount = match hk_count_locks.4 {
        LockPosition::ThisLock => hk_count_locks.2.as_ref().unwrap(),
        LockPosition::CurrentLock => &hk_count_locks.0,
        LockPosition::OtherLock => hk_count_locks.1.as_ref().unwrap(),
    };

    if current_sk.is_canonical_in_the_read() == next_sk.is_canonical_in_the_read() {
        next_hk_count
            .get_extended_hyperkmer_left_id(
                hyperkmers,
                &large_hyperkmers,
                &next_minimizer,
                right_extended_hk,
            )
            .expect("Hash collision on superkmers. Please change your seed.")
    } else {
        next_hk_count
            .get_extended_hyperkmer_right_id(
                hyperkmers,
                &large_hyperkmers,
                &next_minimizer,
                right_extended_hk,
            )
            .expect("Hash collision on superkmers. Please change your seed.")
    }
}

fn is_previous_sk_solid(
    sk_count_locks: &(
        RwLockWriteGuard<SuperKmerCounts>,
        Option<RwLockReadGuard<SuperKmerCounts>>,
        Option<RwLockReadGuard<SuperKmerCounts>>,
        LockPosition,
        LockPosition,
    ),
    previous_sk: &Superkmer,
    threshold: Count,
) -> bool {
    let previous_sk_count: &SuperKmerCounts = match sk_count_locks.3 {
        LockPosition::CurrentLock => &sk_count_locks.0,
        LockPosition::ThisLock => sk_count_locks.1.as_ref().unwrap(),
        LockPosition::OtherLock => sk_count_locks.2.as_ref().unwrap(),
    };
    previous_sk_count.get_count_superkmer(previous_sk) >= threshold
}

fn is_next_sk_solid(
    sk_count_locks: &(
        RwLockWriteGuard<SuperKmerCounts>,
        Option<RwLockReadGuard<SuperKmerCounts>>,
        Option<RwLockReadGuard<SuperKmerCounts>>,
        LockPosition,
        LockPosition,
    ),
    next_sk: &Superkmer,
    threshold: Count,
) -> bool {
    let next_sk_count: &SuperKmerCounts = match sk_count_locks.4 {
        LockPosition::CurrentLock => &sk_count_locks.0,
        LockPosition::ThisLock => sk_count_locks.2.as_ref().unwrap(),
        LockPosition::OtherLock => sk_count_locks.1.as_ref().unwrap(),
    };
    next_sk_count.get_count_superkmer(next_sk) >= threshold
}
// TODO "style" find a better name for the first stage function
#[allow(clippy::too_many_arguments)]
fn first_stage_for_a_chunck(
    sequences: &mut Vec<Vec<u8>>,
    k: usize,
    m: usize,
    threshold: Count,
    sk_count: &Buckets<SuperKmerCounts>,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
) {
    replace_n(sequences);
    let hyperkmers = hyperkmers.read().unwrap();
    for sequence in sequences {
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

            // get sk_count locks
            let sk_count_locks = sk_count.acquire_write_locks(
                current_sk.get_minimizer(),
                previous_sk.get_minimizer(),
                next_sk.get_minimizer(),
            );
            // get hk_count locks
            let hk_count_locks = hk_count.acquire_write_locks(
                current_sk.get_minimizer(),
                previous_sk.get_minimizer(),
                next_sk.get_minimizer(),
            );

            let previous_sk_is_solid =
                is_previous_sk_solid(&sk_count_locks, &previous_sk, threshold);
            let next_sk_is_solid = is_next_sk_solid(&sk_count_locks, &next_sk, threshold);

            // OPTIMIZE est-il possible de stocker les counts pour ne pas les recalculer ?
            let mut current_sk_sount = sk_count_locks.0;
            let current_count = current_sk_sount.increase_count_superkmer(&current_sk);
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
                let (id_left_bucket, id_left_hk, is_large_left) = if previous_sk_is_solid {
                    get_bucket_of_previous_hk(
                        &current_sk,
                        &previous_sk,
                        &hk_count_locks,
                        &left_extended_hk.0,
                        &hyperkmers,
                        large_hyperkmers,
                    )
                } else {
                    // previous sk is not solid => our hyperkmer is not already present
                    // let's add it
                    add_new_hyperkmer(
                        left_extended_hk.3,
                        &left_extended_hk.0,
                        &hyperkmers,
                        large_hyperkmers,
                    )
                };
                debug_assert!(id_left_bucket < NB_BUCKETS);
                let (id_right_bucket, id_right_hk, is_large_right) = if next_sk_is_solid {
                    get_bucket_of_next_hk(
                        &current_sk,
                        &next_sk,
                        &hk_count_locks,
                        &right_extended_hk.0,
                        &hyperkmers,
                        large_hyperkmers,
                    )
                } else {
                    // previous sk is not solid => our hyperkmer is not already present
                    // let's add it
                    add_new_hyperkmer(
                        right_extended_hk.3,
                        &right_extended_hk.0,
                        &hyperkmers,
                        large_hyperkmers,
                    )
                };
                debug_assert!(id_right_bucket < NB_BUCKETS);

                // TODO drop sk before ?
                drop(sk_count_locks.1);
                drop(sk_count_locks.2);
                drop(hk_count_locks.1);
                drop(hk_count_locks.2);

                // we have two ids (left and rigth) of extended hyperkmers containing our left and right hyperkemr
                // let's get their orientation wrt to the orientation of the minimizer

                let left_change_orientation = !left_extended_hk.0.is_canonical();
                let right_change_orientation = !right_extended_hk.0.is_canonical();

                let left_hk_metadata = HKMetadata::new(
                    id_left_bucket,
                    id_left_hk,
                    left_extended_hk.1,
                    left_extended_hk.2,
                    is_large_left,
                    left_change_orientation,
                );

                let right_hk_metadata = HKMetadata::new(
                    id_right_bucket,
                    id_right_hk,
                    right_extended_hk.1,
                    right_extended_hk.2,
                    is_large_right,
                    right_change_orientation,
                );

                #[cfg(debug_assertions)]
                check_correct_inclusion_first_stage(
                    m,
                    &hyperkmers,
                    large_hyperkmers,
                    left_change_orientation,
                    &left_hk_metadata,
                    &left_extended_hk,
                    right_change_orientation,
                    &right_hk_metadata,
                    &right_extended_hk,
                );

                let mut current_hk_count = hk_count_locks.0;
                // TODO
                // let current_hk_count2 = hk_count_locks[cur_hk_count_lock_pos].unwrap();
                current_hk_count.insert_new_entry_in_hyperkmer_count(
                    &current_sk.get_minimizer(),
                    &left_hk_metadata,
                    &right_hk_metadata,
                    current_count,
                );
            } else if current_count > threshold {
                // TODO fusionnner les deux passes
                let large_hyperkmers = large_hyperkmers.read().unwrap();
                let (left_sk, right_sk) = get_left_and_rigth_of_sk(&current_sk);
                let mut current_hk_count = hk_count_locks.0;
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
}

// TODO "style" find a better name for the second stage function
#[allow(clippy::too_many_arguments)]
fn second_stage_for_a_chunk(
    sk_count: &Buckets<SuperKmerCounts>,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &Arc<RwLock<ParallelExtendedHyperkmers>>,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
    sequences: &mut Vec<Vec<u8>>, // OPTIMIZE prendre un iterateur sur des &[u8] ?
    k: usize,
    m: usize,
    threshold: Count,
    discarded_minimizers: &Buckets<HashMap<Minimizer, Count>>,
) {
    replace_n(sequences);
    // large hyperkmers are not going to be modified
    // => we can acquire them here
    let hyperkmers = hyperkmers.read().unwrap();
    let large_hyperkmers = large_hyperkmers.read().unwrap();
    for sequence in sequences {
        let mut superkmers = match compute_superkmers_linear_streaming(sequence, k, m) {
            Some(superkmers_iter) => superkmers_iter.peekable(), // peekable needed to detect the last superkmer of a read
            None => continue,
        };

        {
            // The first stage skips the first superkmer.
            // Therefore, we have to check if the superkmer is present in the index (it might have been indexed if it is in the middle of a read).
            // If so, increase is count.
            // Otherwise, we insert it with a count of 1.

            // extract the first superkmer, skip if none
            let first_truncated_sk = match superkmers.next() {
                Some(sk) => sk,
                None => continue,
            };

            increase_count_of_sk_or_insert_it(
                k,
                &first_truncated_sk,
                hk_count,
                &hyperkmers,
                &large_hyperkmers,
            );
        }

        let mut last_superkmer = None;
        // iterates over all the superkmers exepct the first one
        while let Some(superkmer) = superkmers.next() {
            if superkmers.peek().is_none() {
                last_superkmer = Some(superkmer);
                break;
            }
            let minimizer = superkmer.get_minimizer();
            let sk_count = sk_count.get_from_id_u64(minimizer);
            let sk_count = sk_count.read().unwrap();

            if sk_count.get_count_superkmer(&superkmer) >= threshold {
                continue;
            }

            let hk_count_lock = hk_count.get_from_id_u64(minimizer);
            let hk_count = hk_count_lock.read().unwrap();

            if !hk_count.contains_minimizer(&minimizer) {
                let discarded_minimizers = discarded_minimizers.get_from_id_u64(minimizer);
                let mut discarded_minimizers = discarded_minimizers.write().unwrap();
                // increase count, set to 1 if it was 0
                let count = discarded_minimizers
                    .entry(minimizer)
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

            drop(hk_count);
            let mut hk_count = hk_count_lock.write().unwrap();
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

        {
            // The first stage skips the last superkmer.
            // Therefore, we have to check if the superkmer is present in the index (it might have been indexed if it is in the middle of a read).
            // If so, increase is count.
            // Otherwise, we insert it with a count of 1.
            if let Some(last_truncated_superkmer) = last_superkmer {
                increase_count_of_sk_or_insert_it(
                    k,
                    &last_truncated_superkmer,
                    hk_count,
                    &hyperkmers,
                    &large_hyperkmers,
                );
            }
        }
    }
}

fn increase_count_of_sk_or_insert_it(
    k: usize,
    superkmer: &Superkmer,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &[LargeExtendedHyperkmer],
) {
    let minimizer = superkmer.get_minimizer();

    // get locks
    let hk_count_lock = hk_count.get_from_id_u64(minimizer);
    let mut hk_count = hk_count_lock.write().unwrap();

    // get left and right parts of the superkmer
    let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

    hk_count.increase_count_of_sk_if_found_else_insert_it(
        k,
        hyperkmers,
        large_hyperkmers,
        &minimizer,
        &left_sk,
        &right_sk,
    );
    #[cfg(debug_assertions)]
    {
        debug_assert!(hk_count.search_exact_match(
            &minimizer,
            hyperkmers,
            large_hyperkmers,
            &left_sk,
            &right_sk
        ))
    }
}

#[cfg(debug_assertions)]
fn check_correct_inclusion_first_stage(
    m: usize,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
    left_change_orientation: bool,
    left_hk_metadata: &HKMetadata,
    left_extended_hk: &(Subsequence<NoBitPacked<'_>>, usize, usize, bool),
    right_change_orientation: bool,
    right_hk_metadata: &HKMetadata,
    right_extended_hk: &(Subsequence<NoBitPacked<'_>>, usize, usize, bool),
) {
    use crate::index::components::get_subsequence_from_metadata;

    let (left_hyperkmers, right_hyperkmers) = hyperkmers.acquire_two_locks_read_mode(
        left_hk_metadata.get_bucket_id(),
        right_hk_metadata.get_bucket_id(),
    );
    let right_hyperkmers = match right_hyperkmers.as_ref() {
        Some(x) => x,
        None => &left_hyperkmers,
    };

    let large_hyperkmers = large_hyperkmers.read().unwrap();

    let candidate_left_ext_hk =
        &get_subsequence_from_metadata(&left_hyperkmers, &large_hyperkmers, left_hk_metadata)
            .change_orientation_if(left_change_orientation);
    let candidate_right_ext_hk =
        &get_subsequence_from_metadata(right_hyperkmers, &large_hyperkmers, right_hk_metadata)
            .change_orientation_if(right_change_orientation);

    debug_assert!(left_extended_hk.0.equal_bitpacked(candidate_left_ext_hk));
    debug_assert!(right_extended_hk.0.equal_bitpacked(candidate_right_ext_hk));
    debug_assert!(left_extended_hk.0.to_canonical().equal_bitpacked(
        &get_subsequence_from_metadata(&left_hyperkmers, &large_hyperkmers, left_hk_metadata,)
    ));
    debug_assert!(right_extended_hk.0.to_canonical().equal_bitpacked(
        &get_subsequence_from_metadata(right_hyperkmers, &large_hyperkmers, right_hk_metadata,)
    ));

    let left_ext_hk =
        get_subsequence_from_metadata(&left_hyperkmers, &large_hyperkmers, left_hk_metadata)
            .change_orientation_if(left_hk_metadata.get_change_orientation());

    let right_ext_hk =
        get_subsequence_from_metadata(right_hyperkmers, &large_hyperkmers, right_hk_metadata)
            .change_orientation_if(right_hk_metadata.get_change_orientation());

    // extract candidate hyperkmers
    let left_hyperkmer =
        &left_ext_hk.subsequence(left_hk_metadata.get_start(), left_hk_metadata.get_end());
    let right_hyperkmer =
        &right_ext_hk.subsequence(right_hk_metadata.get_start(), right_hk_metadata.get_end());

    let left_string = left_hyperkmer.to_string();
    let right_string = right_hyperkmer.to_string();

    debug_assert_eq!(
        left_string[(left_string.len() - (m - 2))..left_string.len()],
        right_string[0..(m - 2)]
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_replace_n() {
        let input0 = String::from("ATGGAGCAGCTGACGANNNATGCA");
        let input1 = String::from("ANNAGTAGNCNGAT");
        let input2 = String::from("NNACNGATTAGCN");

        let expected = vec![
            String::from("ATGGAGCAGCTGACGAAAAATGCA"),
            String::from("AAAAGTAGACAGAT"),
            String::from("AAACAGATTAGCA"),
        ];
        let mut inputs = vec![
            input0.as_bytes().to_vec(),
            input1.as_bytes().to_vec(),
            input2.as_bytes().to_vec(),
        ];

        replace_n(&mut inputs);

        let got = inputs
            .into_iter()
            .map(|v| String::from_utf8(v).unwrap())
            .collect_vec();

        assert_eq!(expected, got);
    }
}
