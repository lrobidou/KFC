use crate::get_start_of_minimizer_in_superkmer;
use crate::reverse_complement;
// use crate::PosOfHyperkmerInExtHyperkmer;
use crate::SuperKmerInfos;

pub fn get_left_and_rigth_hk(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> (String, String) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

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
            // we work on the rc of the superkmer
            // as distances are given as if the superkmer is canonical in the read, we must do some math
            assert!(current_sk.start_of_minimizer_as_read < previous_sk.start_of_minimizer_as_read);
            assert!(current_sk.start_of_minimizer_as_read > next_sk.start_of_minimizer_as_read);

            let start_of_minimizer_in_sk = current_sk.superkmer.len()
                - current_sk.start_of_minimizer_as_read
                + current_sk.start_of_superkmer_as_read
                - m;
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
    let end_of_left_hk = start_of_minimizer_in_sk + m - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;
    let end_of_right_hk = start_of_minimizer_in_sk + distance_to_right + m - 1;

    // this is the sequence in between minimizers, including m-1 bases of the minimizers
    let left_hk = &current_sk.superkmer[start_of_left_hk..end_of_left_hk];
    let right_hk = &current_sk.superkmer[start_of_right_hk..end_of_right_hk];

    (String::from(left_hk), String::from(right_hk))
}

/// return the left and right extended hyperkmers of the currrent superkmer
/// extended left hyperkmers is the largest sequence between:
/// - the left part of a superkmer
/// - the right part of the left superkmer
/// returned extended hyerkmer are not in canonical form, but are as they appear in the current superkmer
/// (i.e. as if the the current superkmer was in its canonical form in the read)
/// returned tuple is (left, right)
pub fn get_left_and_rigth_extended_hk(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> ((String, usize, usize), (String, usize, usize)) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

    let start_of_minimizer_in_sk = get_start_of_minimizer_in_superkmer(current_sk);
    let end_of_left_hk = start_of_minimizer_in_sk + m - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;

    // this is the left and rigth part of the superkmer
    let current_left_sk = &current_sk.superkmer[0..end_of_left_hk];
    let current_right_sk = &current_sk.superkmer[start_of_right_hk..current_sk.superkmer.len()];

    let previous_right_sk = if previous_sk.was_read_canonical == current_sk.was_read_canonical {
        String::from(
            &previous_sk.superkmer
                [get_start_of_minimizer_in_superkmer(previous_sk) + 1..previous_sk.superkmer.len()],
        )
    } else {
        reverse_complement(
            &previous_sk.superkmer[0..get_start_of_minimizer_in_superkmer(previous_sk) + m - 1],
        )
    };
    let next_left_sk = if next_sk.was_read_canonical == current_sk.was_read_canonical {
        String::from(&next_sk.superkmer[0..get_start_of_minimizer_in_superkmer(next_sk) + m - 1])
    } else {
        reverse_complement(
            &next_sk.superkmer
                [get_start_of_minimizer_in_superkmer(next_sk) + 1..next_sk.superkmer.len()],
        )
    };

    let len_match_left = std::cmp::min(current_left_sk.len(), previous_right_sk.len());
    let len_match_right = std::cmp::min(current_right_sk.len(), next_left_sk.len());

    let extended_left_sk = if current_left_sk.len() > previous_right_sk.len() {
        (
            current_left_sk.into(),
            current_left_sk.len() - len_match_left,
            current_left_sk.len(),
        )
    } else {
        (previous_right_sk.clone(), 0, len_match_left)
    };

    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        (current_right_sk.into(), 0, len_match_right)
    } else {
        (
            next_left_sk.clone(),
            next_left_sk.len() - len_match_right,
            len_match_right,
        )
    };

    (extended_left_sk, extended_right_sk)
}
