use crate::superkmer::{NoBitPacked, SubsequenceMetadata};
use crate::Superkmer;

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
