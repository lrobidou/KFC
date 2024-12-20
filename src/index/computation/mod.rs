use crate::{
    buckets::{Buckets, LockPosition, ThreeLocks, NB_BUCKETS},
    complexity,
    compute_left_and_right::{get_left_and_rigth_extended_hk, get_left_and_rigth_of_sk},
    index::components::{
        search_exact_hyperkmer_match, AllHyperkmerParts, HKCount, HKMetadata, SuperKmerCounts,
    },
    read_modification::replace_n,
    subsequence::{NoBitPacked, Subsequence},
    superkmer::Superkmer,
    superkmers_computation::compute_superkmers_linear_streaming,
    Count, Minimizer,
};

use std::{
    collections::HashMap,
    path::Path,
    sync::{Arc, Mutex},
};

mod cache;
mod lines_iter;

use cache::CachedValue;

use lines_iter::LinesIter;

const COMPLEXITY_THRESHOLD: u16 = 0;

use itertools::Itertools;
use log::warn;
use needletail::parse_fastx_file;

const BATCH_SIZE: usize = 100;

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
    AllHyperkmerParts,
) {
    let sk_count = Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new);
    let hk_count = Buckets::<HKCount>::new(HKCount::new);
    let hyperkmers = AllHyperkmerParts::new(k, 1000);

    let sequences = parse_fastx_file(path).unwrap();
    let sequences = LinesIter::new(sequences);

    rayon::scope(|s| {
        let chunks = sequences.into_iter().chunks(BATCH_SIZE);
        let hyperkmers = Arc::new(&hyperkmers);
        for chunk in chunks.into_iter() {
            let sk_count = sk_count.clone();
            let hk_count = hk_count.clone();
            let lines = chunk.into_iter().collect_vec(); //TODO copy
            let lines = Arc::new(Mutex::new(lines));
            let hyperkmers = hyperkmers.clone();

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
                )
            });
        }
    });
    (sk_count, hk_count, hyperkmers)
}

/// Second stage of the construction of the KFC index.
///
/// Splits each line of the input into chunks, then performs the second stage on each chunk in parallel.  
#[allow(clippy::too_many_arguments)]
pub fn second_stage<P: AsRef<Path>>(
    sk_count: &Buckets<SuperKmerCounts>,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &AllHyperkmerParts,
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
            let discarded_minimizers = discarded_minimizers.clone();
            let lines = chunk.into_iter().collect_vec(); // TODO copy
            let lines = Arc::new(Mutex::new(lines));

            s.spawn(move |_| {
                let mut sequences = lines.lock().unwrap();
                second_stage_for_a_chunk(
                    &sk_count,
                    &hk_count,
                    hyperkmers,
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

fn get_bucket_of_previous_hk(
    current_sk: &Superkmer,
    previous_sk: &Superkmer,
    hk_count_locks: &ThreeLocks<HKCount>,
    left_extended_hk: &Subsequence<NoBitPacked>,
    hyperkmers: &AllHyperkmerParts,
) -> (usize, usize, bool) {
    // used to access buckets and compute hyperkmer id
    // read hk_count table associated with the previous minimizer
    let previous_minimizer = previous_sk.get_minimizer();

    let previous_hk_count: &HKCount = match hk_count_locks.3 {
        LockPosition::ThisLock => hk_count_locks.1.as_ref().unwrap(),
        LockPosition::CurrentLock => &hk_count_locks.0,
        LockPosition::OtherLock => hk_count_locks.2.as_ref().unwrap(),
    };

    if previous_sk.is_canonical_in_the_read() == current_sk.is_canonical_in_the_read() {
        previous_hk_count
            .get_extended_hyperkmer_right_id(hyperkmers, &previous_minimizer, left_extended_hk)
            .expect("Hash collision on superkmers. Please change your seed.")
    } else {
        previous_hk_count
            .get_extended_hyperkmer_left_id(hyperkmers, &previous_minimizer, left_extended_hk)
            .expect("Hash collision on superkmers. Please change your seed.")
    }
}

fn get_bucket_of_next_hk(
    current_sk: &Superkmer,
    next_sk: &Superkmer,
    hk_count_locks: &ThreeLocks<HKCount>,
    right_extended_hk: &Subsequence<NoBitPacked>,
    hyperkmers: &AllHyperkmerParts,
) -> (usize, usize, bool) {
    let next_minimizer = next_sk.get_minimizer();
    let next_hk_count: &HKCount = match hk_count_locks.4 {
        LockPosition::ThisLock => hk_count_locks.2.as_ref().unwrap(),
        LockPosition::CurrentLock => &hk_count_locks.0,
        LockPosition::OtherLock => hk_count_locks.1.as_ref().unwrap(),
    };

    if current_sk.is_canonical_in_the_read() == next_sk.is_canonical_in_the_read() {
        next_hk_count
            .get_extended_hyperkmer_left_id(hyperkmers, &next_minimizer, right_extended_hk)
            .expect("Hash collision on superkmers. Please change your seed.")
    } else {
        next_hk_count
            .get_extended_hyperkmer_right_id(hyperkmers, &next_minimizer, right_extended_hk)
            .expect("Hash collision on superkmers. Please change your seed.")
    }
}

fn is_previous_sk_solid(
    sk_count_locks: &ThreeLocks<SuperKmerCounts>,
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
    sk_count_locks: &ThreeLocks<SuperKmerCounts>,
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

fn get_bucket_of_previous_hk_or_insert_if_not_found(
    current_sk: &Superkmer,
    previous_sk: &Superkmer,
    previous_sk_is_solid: bool,
    hk_count_locks: &ThreeLocks<HKCount>,
    hyperkmers: &AllHyperkmerParts,
    left_extended_hk: (Subsequence<NoBitPacked>, usize, usize, bool),
    cached_value: &Option<CachedValue>,
) -> (usize, usize, bool) {
    if previous_sk_is_solid {
        if current_sk.is_canonical_in_the_read() {
            debug_assert!(previous_sk.start_of_minimizer() < current_sk.start_of_minimizer());
            // previous sk might be inserted already, let's look at the cache
            if let Some(cached_value) = cached_value {
                return (
                    cached_value.get_id_bucket(),
                    cached_value.get_id_hk(),
                    cached_value.get_is_large(),
                );
                // let u = get_bucket_of_previous_hk(
                //     current_sk,
                //     previous_sk,
                //     hk_count_locks,
                //     &left_extended_hk.0,
                //     hyperkmers,
                //     large_hyperkmers,
                // );
                // assert_eq!(v, u);
            }
        }
        get_bucket_of_previous_hk(
            current_sk,
            previous_sk,
            hk_count_locks,
            &left_extended_hk.0,
            hyperkmers,
        )
    } else {
        // previous sk is not solid => our hyperkmer is not already present
        // let's add it
        hyperkmers.add_new_hyperkmer(left_extended_hk.3, &left_extended_hk.0)
    }
}

fn get_bucket_of_next_hk_or_insert_if_not_found(
    current_sk: &Superkmer,
    next_sk: &Superkmer,
    next_sk_is_solid: bool,
    hk_count_locks: &ThreeLocks<HKCount>,
    hyperkmers: &AllHyperkmerParts,
    right_extended_hk: &(Subsequence<NoBitPacked>, usize, usize, bool),
    cached_value: &Option<CachedValue>,
) -> (usize, usize, bool) {
    if next_sk_is_solid {
        if !current_sk.is_canonical_in_the_read() {
            debug_assert!(next_sk.start_of_minimizer() < current_sk.start_of_minimizer());
            // previous sk might be inserted already, let's look at the cache
            if let Some(cached_value) = cached_value {
                return (
                    cached_value.get_id_bucket(),
                    cached_value.get_id_hk(),
                    cached_value.get_is_large(),
                );
                // let u = get_bucket_of_previous_hk(
                //     current_sk,
                //     previous_sk,
                //     hk_count_locks,
                //     &left_extended_hk.0,
                //     hyperkmers,
                //     large_hyperkmers,
                // );
                // assert_eq!(v, u);
            }
        }

        get_bucket_of_next_hk(
            current_sk,
            next_sk,
            hk_count_locks,
            &right_extended_hk.0,
            hyperkmers,
        )
    } else {
        // previous sk is not solid => our hyperkmer is not already present
        // let's add it
        hyperkmers.add_new_hyperkmer(right_extended_hk.3, &right_extended_hk.0)
    }
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
    hyperkmers: &AllHyperkmerParts,
) {
    replace_n(sequences);
    for sequence in sequences {
        let superkmers = match compute_superkmers_linear_streaming(sequence, k, m) {
            Some(superkmers_iter) => superkmers_iter,
            None => continue, // no superkmer in this string => we skip it
        };
        let mut cached_value: Option<CachedValue> = None;
        for (previous_sk, current_sk, next_sk) in superkmers.into_iter().tuple_windows() {
            // skip superkmer is they only have k-mer with low complexity
            #[allow(clippy::absurd_extreme_comparisons)]
            if COMPLEXITY_THRESHOLD > 0 {
                // only computes complexity if it makes sense
                if !complexity::is_complexity_above_threshold(
                    current_sk.superkmer.get_subsequence_as_in_read(),
                    k,
                    COMPLEXITY_THRESHOLD,
                ) {
                    continue;
                }
            }

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
                let (id_left_bucket, id_left_hk, is_large_left) =
                    get_bucket_of_previous_hk_or_insert_if_not_found(
                        &current_sk,
                        &previous_sk,
                        previous_sk_is_solid,
                        &hk_count_locks,
                        hyperkmers,
                        left_extended_hk,
                        &cached_value,
                    );
                debug_assert!(id_left_bucket < NB_BUCKETS);
                let (id_right_bucket, id_right_hk, is_large_right) =
                    get_bucket_of_next_hk_or_insert_if_not_found(
                        &current_sk,
                        &next_sk,
                        next_sk_is_solid,
                        &hk_count_locks,
                        hyperkmers,
                        &right_extended_hk,
                        &cached_value,
                    );
                debug_assert!(id_right_bucket < NB_BUCKETS);

                // TODO drop sk before ?
                drop(sk_count_locks.1);
                drop(sk_count_locks.2);
                drop(hk_count_locks.1);
                drop(hk_count_locks.2);

                // update the cache
                cached_value = Some(if current_sk.is_canonical_in_the_read() {
                    CachedValue::new(id_right_bucket, id_right_hk, is_large_right)
                } else {
                    CachedValue::new(id_left_bucket, id_left_hk, is_large_left)
                });

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
                    hyperkmers,
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
                // TODO cache
                cached_value = None;
                // TODO fusionnner les deux passes
                let (left_sk, right_sk) = get_left_and_rigth_of_sk(&current_sk);
                let mut current_hk_count = hk_count_locks.0;
                let found = current_hk_count.increase_count_if_exact_match(
                    &current_sk.get_minimizer(),
                    hyperkmers,
                    &left_sk,
                    &right_sk,
                );
                if !found {
                    // if no exact match, then we must at least have an approximate match
                    let new_left_and_right_metadata = current_hk_count.search_for_inclusion(
                        hyperkmers,
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
                        hyperkmers,
                        &left_sk,
                        &right_sk,
                        &metadata_to_insert_left,
                        &metadata_to_insert_right
                    ));
                }
            } else {
                // nothing to report => clear cache
                cached_value = None;
            }
        }
    }
}

// TODO "style" find a better name for the second stage function
#[allow(clippy::too_many_arguments)]
fn second_stage_for_a_chunk(
    sk_count: &Buckets<SuperKmerCounts>,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &AllHyperkmerParts,
    sequences: &mut Vec<Vec<u8>>, // OPTIMIZE prendre un iterateur sur des &[u8] ?
    k: usize,
    m: usize,
    threshold: Count,
    discarded_minimizers: &Buckets<HashMap<Minimizer, Count>>,
) {
    // TODO cache
    replace_n(sequences);
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

            #[allow(clippy::absurd_extreme_comparisons)]
            if COMPLEXITY_THRESHOLD > 0 {
                // only computes complexity if it makes sense
                if !complexity::is_complexity_above_threshold(
                    first_truncated_sk.superkmer.get_subsequence_as_in_read(),
                    k,
                    COMPLEXITY_THRESHOLD,
                ) {
                    continue;
                }
            }

            increase_count_of_sk_or_insert_it(k, &first_truncated_sk, hk_count, hyperkmers);
        }

        let mut last_superkmer = None;
        // iterates over all the superkmers exepct the first one
        while let Some(superkmer) = superkmers.next() {
            if superkmers.peek().is_none() {
                last_superkmer = Some(superkmer);
                break;
            }
            #[allow(clippy::absurd_extreme_comparisons)]
            if COMPLEXITY_THRESHOLD > 0 {
                // only computes complexity if it makes sense
                if !complexity::is_complexity_above_threshold(
                    superkmer.superkmer.get_subsequence_as_in_read(),
                    k,
                    COMPLEXITY_THRESHOLD,
                ) {
                    continue;
                }
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
            let match_metadata = hk_count
                .search_for_maximal_inclusion(hyperkmers, k, m, &minimizer, &left_sk, &right_sk);

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
                #[allow(clippy::absurd_extreme_comparisons)]
                if COMPLEXITY_THRESHOLD > 0 {
                    // only computes complexity if it makes sense
                    if !complexity::is_complexity_above_threshold(
                        last_truncated_superkmer
                            .superkmer
                            .get_subsequence_as_in_read(),
                        k,
                        COMPLEXITY_THRESHOLD,
                    ) {
                        continue;
                    }
                }
                increase_count_of_sk_or_insert_it(
                    k,
                    &last_truncated_superkmer,
                    hk_count,
                    hyperkmers,
                );
            }
        }
    }
}

fn increase_count_of_sk_or_insert_it(
    k: usize,
    superkmer: &Superkmer,
    hk_count: &Buckets<HKCount>,
    hyperkmers: &AllHyperkmerParts,
) {
    let minimizer = superkmer.get_minimizer();

    // get locks
    let hk_count_lock = hk_count.get_from_id_u64(minimizer);
    let mut hk_count = hk_count_lock.write().unwrap();

    // get left and right parts of the superkmer
    let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

    hk_count.increase_count_of_sk_if_found_else_insert_it(
        k, hyperkmers, &minimizer, &left_sk, &right_sk,
    );
    #[cfg(debug_assertions)]
    {
        debug_assert!(hk_count.search_exact_match(&minimizer, hyperkmers, &left_sk, &right_sk))
    }
}

#[cfg(debug_assertions)]
fn check_correct_inclusion_first_stage(
    m: usize,
    hyperkmers: &AllHyperkmerParts,
    left_change_orientation: bool,
    left_hk_metadata: &HKMetadata,
    left_extended_hk: &(Subsequence<NoBitPacked<'_>>, usize, usize, bool),
    right_change_orientation: bool,
    right_hk_metadata: &HKMetadata,
    right_extended_hk: &(Subsequence<NoBitPacked<'_>>, usize, usize, bool),
) {
    let candidate_left_ext_hk = hyperkmers
        .get_subsequence_from_metadata(left_hk_metadata)
        .change_orientation_if(left_change_orientation);
    let candidate_right_ext_hk = hyperkmers
        .get_subsequence_from_metadata(right_hk_metadata)
        .change_orientation_if(right_change_orientation);

    debug_assert!(left_extended_hk.0.equal_bitpacked(&candidate_left_ext_hk));
    debug_assert!(right_extended_hk.0.equal_bitpacked(&candidate_right_ext_hk));
    debug_assert!(left_extended_hk
        .0
        .to_canonical()
        .equal_bitpacked(&hyperkmers.get_subsequence_from_metadata(left_hk_metadata)));
    debug_assert!(right_extended_hk
        .0
        .to_canonical()
        .equal_bitpacked(&hyperkmers.get_subsequence_from_metadata(right_hk_metadata)));

    let left_ext_hk = hyperkmers
        .get_subsequence_from_metadata(left_hk_metadata)
        .change_orientation_if(left_hk_metadata.get_change_orientation());

    let right_ext_hk = hyperkmers
        .get_subsequence_from_metadata(right_hk_metadata)
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
