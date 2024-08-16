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
///     - the left part of a superkmer
///     - the right part of the left superkmer
/// returned extended hyerkmer are not in canonical form, but are as they appear in the current superkmer
/// (i.e. as if the the current superkmer was in its canonical form in the read)
/// returned tuple is (left, right)
pub fn get_left_and_rigth_extended_hk<'a>(
    previous_sk: &Superkmer<'a>,
    current_sk: &Superkmer<'a>,
    next_sk: &Superkmer<'a>,
    k: usize,
) -> (
    (SubsequenceMetadata<'a, NoBitPacked>, usize, usize, bool),
    (SubsequenceMetadata<'a, NoBitPacked>, usize, usize, bool),
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

    let extended_left_sk = if current_left_sk.len() > previous_right_sk.len() {
        (current_left_sk, 0, current_left_sk.len(), false)
    } else {
        (previous_right_sk, 0, current_left_sk.len(), false)
    };

    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        (current_right_sk, 0, current_right_sk.len(), false)
    } else {
        (
            next_left_sk,
            next_left_sk.len() - current_right_sk.len(),
            next_left_sk.len(),
            false,
        )
    };

    debug_assert!(extended_left_sk.0.len() < k);
    debug_assert!(extended_right_sk.0.len() < k);

    // TODO discuss should the test be about if the two minimizers are the same?
    // TODO document this code
    // I (lucas) feel like I am missing something here
    // I think there's an edge case and the "inclusion hypothesis" can be violated if the two minimizers are equal
    let extended_left_sk = if extended_left_sk.0.len() == k - 1 {
        extended_left_sk
    } else if unlikely(extended_left_sk.0.len() < k - 1) {
        // the right context is not maximal => need to create such context
        if current_sk.is_canonical_in_the_read() {
            // TODO review my code please
            debug_assert!(previous_sk.superkmer.start() < current_sk.superkmer.start());

            let left_ext = SubsequenceMetadata::new(
                previous_sk.read,
                previous_sk.superkmer.start(),
                current_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );

            let start = left_ext.len() - current_sk.superkmer.len();
            (left_ext, start, start + current_left_sk.len(), true)
        } else {
            // TODO review my code please
            debug_assert!(previous_sk.superkmer.start() > current_sk.superkmer.start());

            let left_ext = SubsequenceMetadata::new(
                current_sk.read,
                current_sk.superkmer.start(),
                previous_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );

            let start = left_ext.len() - current_sk.superkmer.len();
            (left_ext, start, start + current_left_sk.len(), true)
        }
    } else {
        unreachable!(
            "the left context of the superkmer {} is larger or equal to k",
            current_sk.superkmer
        )
    };

    #[cfg(debug_assertions)]
    {
        // left hyperkmer should end by the minimizer (excluding its last character)
        let left_hyperkmer =
            extended_left_sk.0.to_string()[extended_left_sk.1..extended_left_sk.2].to_string(); // potentially its revcomp
        let m = current_sk.end_of_minimizer() - current_sk.start_of_minimizer();
        let minimizer = current_sk.minimizer_string();
        debug_assert_eq!(
            left_hyperkmer[left_hyperkmer.len() + 1 - m..left_hyperkmer.len()].to_owned(),
            minimizer[0..minimizer.len() - 1]
        );
    }

    // TODO discuss (same as above)
    let extended_right_sk = if extended_right_sk.0.len() == k - 1 {
        extended_right_sk
    } else if unlikely(extended_right_sk.0.len() < k - 1) {
        // the right context is not maximal => need to create such context
        if current_sk.is_canonical_in_the_read() {
            // TODO review my code please
            debug_assert!(current_sk.superkmer.start() < (next_sk.superkmer.start()));

            let right_ext = SubsequenceMetadata::new(
                current_sk.read,
                current_sk.superkmer.start(),
                next_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );

            let m = current_sk.end_of_minimizer() - current_sk.start_of_minimizer();
            let start = current_left_sk.len() - m + 2;
            (right_ext, start, start + current_right_sk.len(), true)
        } else {
            // TODO review my code please
            debug_assert!(current_sk.superkmer.start() > (next_sk.superkmer.start()));

            let right_ext = SubsequenceMetadata::new(
                next_sk.read,
                next_sk.superkmer.start(),
                current_sk.superkmer.end(),
                current_sk.is_canonical_in_the_read(),
            );

            let m = current_sk.end_of_minimizer() - current_sk.start_of_minimizer();
            let start = current_left_sk.len() - m + 2;
            (right_ext, start, start + current_right_sk.len(), true)
        }
    } else {
        unreachable!(
            "the right context of the superkmer {} is larger or equal to k",
            current_sk.superkmer
        )
    };

    #[cfg(debug_assertions)]
    {
        // right hyperkmer should start by the minimizer (excluding its first character)
        let right_hyperkmer =
            extended_right_sk.0.to_string()[extended_right_sk.1..extended_right_sk.2].to_string();
        let m = current_sk.end_of_minimizer() - current_sk.start_of_minimizer();
        let minimizer = current_sk.minimizer_string();
        debug_assert_eq!(right_hyperkmer[0..m - 1].to_owned(), minimizer[1..m]);
    }
    if !extended_left_sk.3 {
        debug_assert!(extended_left_sk.0.len() == k - 1);
    }
    if !extended_right_sk.3 {
        debug_assert!(extended_right_sk.0.len() == k - 1);
    }
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
