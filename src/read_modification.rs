//! Modifies the input sequences by removing `N`s or homopolymers

// TODO find better name for this module

pub fn replace_n(sequences: &mut Vec<Vec<u8>>) {
    for sequence in sequences {
        for char in sequence {
            if char == &b'N' {
                *char = b'A';
            }
        }
    }
}

// TODO test
pub fn remove_homopolymers_and_n(sequences: &mut Vec<Vec<u8>>) {
    for sequence in sequences {
        let n = sequence.len();
        let mut i = 1;
        let mut last_chr = if sequence[0] == b'N' {
            b'A'
        } else {
            sequence[0]
        };
        sequence[0] = last_chr;
        for pos_in_original in 1..n {
            let current_chr = sequence[pos_in_original];
            let current_chr = if current_chr == b'N' {
                b'A'
            } else {
                current_chr
            };
            if last_chr != current_chr {
                last_chr = current_chr;
                sequence[i] = current_chr;
                i += 1;
            }
        }
        // remove element after the end of the vector
        sequence.truncate(i);
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

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

    #[test]
    fn test_remove_homopolymers_and_n() {
        let input0 = String::from("ATGGAGCAGCTGACGANNNATGCA");
        let input1 = String::from("ANNAGTAGNCNGAT");
        let input2 = String::from("NNACNGATTAGCN");

        let expected = vec![
            String::from("ATGAGCAGCTGACGATGCA"),
            String::from("AGTAGACAGAT"),
            String::from("ACAGATAGCA"),
        ];
        let mut inputs = vec![
            input0.as_bytes().to_vec(),
            input1.as_bytes().to_vec(),
            input2.as_bytes().to_vec(),
        ];

        remove_homopolymers_and_n(&mut inputs);

        let got = inputs
            .into_iter()
            .map(|v| String::from_utf8(v).unwrap())
            .collect_vec();

        assert_eq!(expected, got);
    }
}
