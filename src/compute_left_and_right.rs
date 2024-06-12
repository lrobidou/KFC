use crate::reverse_complement;
use crate::PosOfHyperkmerInExtHyperkmer;
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
) -> (
    (String, PosOfHyperkmerInExtHyperkmer, usize),
    (String, PosOfHyperkmerInExtHyperkmer, usize),
) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

    // get the starting position of the minimizer in the superkmer
    let get_start_of_minimizer = |sk: &SuperKmerInfos| {
        if sk.was_read_canonical {
            sk.start_of_minimizer_as_read - sk.start_of_superkmer_as_read
        } else {
            sk.superkmer.len() + sk.start_of_superkmer_as_read - sk.start_of_minimizer_as_read - m
        }
    };

    let start_of_minimizer_in_sk = get_start_of_minimizer(current_sk);

    let end_of_left_hk = start_of_minimizer_in_sk + m - 1;
    let start_of_right_hk = start_of_minimizer_in_sk + 1;

    // this is the left and rigth part of the superkmer
    let current_left_sk = &current_sk.superkmer[0..end_of_left_hk];
    let current_right_sk = &current_sk.superkmer[start_of_right_hk..current_sk.superkmer.len()];

    let previous_right_sk = if previous_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(
            &previous_sk.superkmer
                [get_start_of_minimizer(previous_sk) + 1..previous_sk.superkmer.len()],
        )
    } else {
        reverse_complement(&previous_sk.superkmer[0..get_start_of_minimizer(previous_sk) + m - 1])
    };
    let next_left_sk = if next_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(&next_sk.superkmer[0..get_start_of_minimizer(next_sk) + m - 1])
    } else {
        reverse_complement(
            &next_sk.superkmer[get_start_of_minimizer(next_sk) + 1..next_sk.superkmer.len()],
        )
    };

    // TODO remove checks
    let extended_left_sk = if current_left_sk.len() > previous_right_sk.len() {
        assert!(current_left_sk.contains(&previous_right_sk));
        assert!(current_left_sk.ends_with(&previous_right_sk));
        (current_left_sk.into(), PosOfHyperkmerInExtHyperkmer::End)
    } else {
        assert!(previous_right_sk.contains(current_left_sk));
        assert!(previous_right_sk.starts_with(current_left_sk));
        (
            previous_right_sk.clone(),
            PosOfHyperkmerInExtHyperkmer::Start,
        )
    };
    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        assert!(current_right_sk.contains(&next_left_sk));
        assert!(current_right_sk.starts_with(&next_left_sk));
        (current_right_sk.into(), PosOfHyperkmerInExtHyperkmer::Start)
    } else {
        assert!(next_left_sk.contains(current_right_sk));
        assert!(next_left_sk.ends_with(current_right_sk));
        (next_left_sk.clone(), PosOfHyperkmerInExtHyperkmer::End)
    };
    (
        (
            extended_left_sk.0,
            extended_left_sk.1,
            std::cmp::min(current_left_sk.len(), previous_right_sk.len()),
        ),
        (
            extended_right_sk.0,
            extended_right_sk.1,
            std::cmp::min(current_right_sk.len(), next_left_sk.len()),
        ),
    )
}

/// same, but with assert
pub fn get_left_and_rigth_extended_hk_dbg(
    previous_sk: &SuperKmerInfos,
    current_sk: &SuperKmerInfos,
    next_sk: &SuperKmerInfos,
) -> (
    (String, PosOfHyperkmerInExtHyperkmer, usize),
    (String, PosOfHyperkmerInExtHyperkmer, usize),
) {
    // Caution: the next and previous superkmer are given as they appear in the read.
    // * but still in the order they would appear if the current superkmer was canonical *
    // this leads to conceptually having to reverse the left and right sequences' content
    // if the superkmer was not read in its canonical form
    let m = current_sk.minimizer.len();

    // get the starting position of the minimizer in the superkmer
    let get_start_of_minimizer = |sk: &SuperKmerInfos| {
        if sk.was_read_canonical {
            sk.start_of_minimizer_as_read - sk.start_of_superkmer_as_read
        } else {
            sk.superkmer.len() + sk.start_of_superkmer_as_read - sk.start_of_minimizer_as_read - m
        }
    };

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
                + current_sk.start_of_superkmer_as_read
                - current_sk.start_of_minimizer_as_read
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

    // this is the left and rigth part of the superkmer
    let current_left_sk = &current_sk.superkmer[0..end_of_left_hk];
    let current_right_sk = &current_sk.superkmer[start_of_right_hk..current_sk.superkmer.len()];

    let previous_right_sk = if previous_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(
            &previous_sk.superkmer
                [get_start_of_minimizer(previous_sk) + 1..previous_sk.superkmer.len()],
        )
    } else {
        reverse_complement(&previous_sk.superkmer[0..get_start_of_minimizer(previous_sk) + m - 1])
    };
    let next_left_sk = if next_sk.was_read_canonical == current_sk.was_read_canonical {
        // reverse left and right
        String::from(&next_sk.superkmer[0..get_start_of_minimizer(next_sk) + m - 1])
    } else {
        reverse_complement(
            &next_sk.superkmer[get_start_of_minimizer(next_sk) + 1..next_sk.superkmer.len()],
        )
    };

    let left_hk = &current_sk.superkmer[start_of_left_hk..end_of_left_hk];
    let right_hk = &current_sk.superkmer[start_of_right_hk..end_of_right_hk];
    // TODO remove checks
    let extended_left_sk = if current_left_sk.len() > previous_right_sk.len() {
        assert!(current_left_sk.contains(&previous_right_sk));
        assert!(current_left_sk.ends_with(&previous_right_sk));
        assert!(current_left_sk.contains(left_hk));
        assert!(current_left_sk.ends_with(left_hk));
        // if current_left_sk == "CTGCCTGATGGAGGGGGATAACTACTGGAA" {
        //     println!("match 1");
        //     println!("previous_sk = {:?}", previous_sk);
        //     println!("current_sk = {:?}", current_sk);
        //     println!("next_sk = {:?}", next_sk);
        //     println!("left_hk = {}", left_hk);
        //     panic!();
        //     // match 1
        //     // left_hk = TCCGGTAACGGACCGAGTTCAGAAA
        // }
        (current_left_sk.into(), PosOfHyperkmerInExtHyperkmer::End)
    } else {
        assert!(previous_right_sk.contains(current_left_sk));
        assert!(previous_right_sk.starts_with(current_left_sk));
        assert!(previous_right_sk.contains(left_hk));
        assert!(previous_right_sk.starts_with(left_hk));
        (
            previous_right_sk.clone(),
            PosOfHyperkmerInExtHyperkmer::Start,
        )
    };
    let extended_right_sk = if current_right_sk.len() > next_left_sk.len() {
        assert!(current_right_sk.contains(&next_left_sk));
        assert!(current_right_sk.starts_with(&next_left_sk));
        assert!(current_right_sk.contains(right_hk));
        assert!(current_right_sk.starts_with(right_hk));
        (current_right_sk.into(), PosOfHyperkmerInExtHyperkmer::Start)
    } else {
        assert!(next_left_sk.contains(current_right_sk));
        assert!(next_left_sk.ends_with(current_right_sk));
        assert!(next_left_sk.contains(right_hk));
        assert!(next_left_sk.ends_with(right_hk));
        if next_left_sk == "ACGGACCGAGTTCAGAAATAAATAACGCGT" {
            println!("match 2");
            println!("right_hk = {}", right_hk);
            // match 2
            // right_hk = ACGGACCGAGTTCAGAAATAAATAACGCGT
        }
        (next_left_sk.clone(), PosOfHyperkmerInExtHyperkmer::End)
    };
    (
        (
            extended_left_sk.0,
            extended_left_sk.1,
            std::cmp::min(current_left_sk.len(), previous_right_sk.len()),
        ),
        (
            extended_right_sk.0,
            extended_right_sk.1,
            std::cmp::min(current_right_sk.len(), next_left_sk.len()),
        ),
    )
}
