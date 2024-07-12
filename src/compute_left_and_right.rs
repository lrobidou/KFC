use crate::superkmer::{NoBitPacked, SubsequenceMetadata};
use crate::Superkmer;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

/// return the left and right extended hyperkmers of the currrent superkmer
/// extended left hyperkmer is the largest sequence between:
/// - the left part of a superkmer
/// - the right part of the left superkmer
/// returned extended hyerkmer are not in canonical form, but are as they appear in the current superkmer
/// (i.e. as if the the current superkmer was in its canonical form in the read)
/// returned tuple is (left, right)
pub fn get_left_and_rigth_extended_hk<'a>(
    previous_sk: &Superkmer<'a>,
    current_sk: &Superkmer<'a>,
    next_sk: &Superkmer<'a>,
    k: usize,
) -> (
    (SubsequenceMetadata<'a, NoBitPacked>, usize, usize),
    (SubsequenceMetadata<'a, NoBitPacked>, usize, usize),
) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    // let m = current_sk.minimizer.len();

    // this is the left and rigth part of the superkmer
    let (previous_left_sk, previous_right_sk) = get_left_and_rigth_of_sk(previous_sk);
    let (current_left_sk, current_right_sk) = get_left_and_rigth_of_sk(current_sk);
    let (next_left_sk, next_right_sk) = get_left_and_rigth_of_sk(next_sk);

    let previous_right_sk =
        if previous_sk.is_canonical_in_the_read() == current_sk.is_canonical_in_the_read() {
            previous_right_sk
        } else {
            previous_left_sk.change_orientation()
        };
    let next_left_sk =
        if next_sk.is_canonical_in_the_read() == current_sk.is_canonical_in_the_read() {
            next_left_sk
        } else {
            next_right_sk.change_orientation()
        };

    let extended_left_sk: (SubsequenceMetadata<'a, NoBitPacked>, usize, usize) =
        if current_left_sk.len() > previous_right_sk.len() {
            (current_left_sk, 0, current_left_sk.len())
        } else {
            (previous_right_sk, 0, current_left_sk.len())
        };

    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        (current_right_sk, 0, current_right_sk.len())
    } else {
        (
            next_left_sk,
            next_left_sk.len() - current_right_sk.len(),
            next_left_sk.len(),
        )
    };

    debug_assert!(extended_left_sk.0.len() < k);
    debug_assert!(extended_right_sk.0.len() < k);

    let extended_left_sk = if extended_left_sk.0.len() == k - 1 {
        extended_left_sk
    } else if unlikely(extended_left_sk.0.len() < k - 1) {
        // the right context is not maximal => need to create such context
        if current_sk.is_canonical_in_the_read() {
            debug_assert!(previous_sk.superkmer.len() >= k);
            let left_ext = SubsequenceMetadata::new(
                previous_sk.read,
                previous_sk.superkmer.end() - (k - 1),
                previous_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );
            (
                left_ext,
                left_ext.len() - current_left_sk.len(),
                left_ext.len(),
            )
        } else {
            debug_assert!(previous_sk.superkmer.len() >= k);
            let left_ext = SubsequenceMetadata::new(
                previous_sk.read,
                previous_sk.superkmer.start(),
                previous_sk.superkmer.start() + (k - 1),
                current_sk.is_canonical_in_the_read(),
            );
            (left_ext, 0, current_left_sk.len())
        }
    } else {
        panic!() // TODO error message
    };

    let extended_right_sk = if extended_right_sk.0.len() == k - 1 {
        extended_right_sk
    } else if unlikely(extended_right_sk.0.len() < k - 1) {
        // the right context is not maximal => need to create such context
        if current_sk.is_canonical_in_the_read() {
            let right_ext = SubsequenceMetadata::new(
                current_sk.read,
                current_sk.superkmer.end() - (k - 1),
                current_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );
            (
                right_ext,
                right_ext.len() - current_right_sk.len(),
                right_ext.len(),
            )
        } else {
            let right_ext = SubsequenceMetadata::new(
                current_sk.read,
                current_sk.superkmer.start(),
                current_sk.superkmer.start() + (k - 1),
                current_sk.is_canonical_in_the_read(),
            );
            (right_ext, 0, current_right_sk.len())
        }
    } else {
        panic!() // TODO error message
    };
    debug_assert!(extended_left_sk.0.len() == k - 1);
    debug_assert!(extended_right_sk.0.len() == k - 1);

    (extended_left_sk, extended_right_sk)
}

// left and right part of the canonical superkmer
pub fn get_left_and_rigth_of_sk<'a>(
    superkmer: &Superkmer<'a>,
) -> (
    SubsequenceMetadata<'a, NoBitPacked>,
    SubsequenceMetadata<'a, NoBitPacked>,
) {
    let left = SubsequenceMetadata::new(
        superkmer.read,
        superkmer.superkmer.start(),
        superkmer.end_of_minimizer() - 1,
        superkmer.is_canonical_in_the_read(),
    );
    let right = SubsequenceMetadata::new(
        superkmer.read,
        superkmer.start_of_minimizer() + 1,
        superkmer.superkmer.end(),
        superkmer.is_canonical_in_the_read(),
    );

    if superkmer.is_canonical_in_the_read() {
        (left, right)
    } else {
        (right, left)
    }
}
