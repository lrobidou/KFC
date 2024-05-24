// Get the reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => base,
        })
        .collect()
}

// Get the canonical form of a k-mer
fn get_canonical_kmer(kmer: &str) -> (bool, String) {
    let rev_comp = reverse_complement(kmer);
    if kmer < &rev_comp {
        (true, kmer.to_string())
    } else {
        (false, rev_comp.to_string())
    }
}

// Get the lexicographically smallest m-mer in a sequence considering canonicity
fn get_canonical_minimizer(seq: &str, m: usize) -> (String, usize, bool) {
    let (mut is_forward, mut minimizer) = get_canonical_kmer(&seq[0..m]);
    let mut position = 0;
    for i in 1..=seq.len() - m {
        let (is_candidate_forward, candidate) = get_canonical_kmer(&seq[i..i + m]);
        if candidate <= minimizer {
            minimizer = candidate;
            is_forward = is_candidate_forward;
            position = i;
        }
    }
    (minimizer, position, is_forward)
}

#[derive(Debug, PartialEq)]
pub struct SuperKmerInfos {
    pub superkmer: String,
    pub minimizer: String,
    pub is_forward: bool,
    pub start_of_minimizer: usize,
    pub start_of_superkmer: usize,
}

// Function to compute superkmers
// OPTIMIZE: this is on O(k*sequence.len()). A solution in O(sequence.len()) is possible
pub fn compute_superkmers(sequence: &str, k: usize, m: usize) -> Vec<SuperKmerInfos> {
    let mut superkmers = Vec::new();
    if sequence.len() < k {
        return superkmers;
    }

    let kmer = &sequence[0..k];
    let (current_minimizer, current_relative_position, forward) = get_canonical_minimizer(kmer, m);
    // superkmer up to here, that will be pushed in the returned vector when the minimizer changes
    let mut superkmerinfos = SuperKmerInfos {
        superkmer: kmer.to_string(), // superkmer will be modified as we see new bases
        minimizer: current_minimizer,
        is_forward: forward,
        start_of_minimizer: current_relative_position,
        start_of_superkmer: 0,
    };

    for i in 1..=sequence.len() - k {
        // for each kmer: if the minimiser the same, append the new base to the superkmer, else create a new superkmer
        let kmer = &sequence[i..i + k];
        let (candidate_minimizer, candidate_relative_minimizer_position, candidate_is_forward) =
            get_canonical_minimizer(kmer, m);
        let candidate_absolute_minimizer_position = candidate_relative_minimizer_position + i;

        // does the minimizer of this kmer have the same value and same position as the minimizer of the previous kmer?
        let is_same_minimizer = candidate_minimizer == superkmerinfos.minimizer;
        let is_at_same_position =
            superkmerinfos.start_of_minimizer == candidate_absolute_minimizer_position;

        if is_same_minimizer && is_at_same_position {
            // if so, add the new base to the superkmer and keep other data untouched
            superkmerinfos
                .superkmer
                .push(sequence.chars().nth(i + k - 1).unwrap());
        } else {
            // otherwise, push the superkmer and create a new one
            superkmers.push(superkmerinfos);
            superkmerinfos = SuperKmerInfos {
                superkmer: kmer.to_string(),
                minimizer: candidate_minimizer,
                is_forward: candidate_is_forward,
                start_of_minimizer: candidate_absolute_minimizer_position,
                start_of_superkmer: i,
            };
        }
    }

    // Add the last superkmer
    superkmers.push(superkmerinfos);

    superkmers
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACTGTGCAGTGCA"), "TGCACTGCACAGT");
    }

    #[test]
    fn test_reverse_complement_n() {
        assert_eq!(reverse_complement("ACTGTGCAGTNNGNCA"), "TGNCNNACTGCACAGT");
    }

    #[test]
    fn test_get_canonical_kmer() {
        assert_eq!(
            get_canonical_kmer("ACTGCGATGACGCAGATAGCAGATAGC"),
            (true, String::from("ACTGCGATGACGCAGATAGCAGATAGC"))
        );
    }

    #[test]
    fn test_get_canonical_kmer2() {
        assert_eq!(
            get_canonical_kmer("GCTATCTGCTATCTGCGTCATCGCAGT"),
            (false, String::from("ACTGCGATGACGCAGATAGCAGATAGC"))
        );
    }

    #[test]
    fn test_get_canonical_minimizer() {
        assert_eq!(
            get_canonical_minimizer("AGCAGCTAGCATTTTTGCAGT", 5),
            (String::from("AAAAA"), 11, false)
        );
    }

    #[test]
    fn test_compute_superkmers() {
        assert_eq!(
            compute_superkmers("AGCAGCTAGCATTTTTGCAGT", 16, 5),
            vec![SuperKmerInfos {
                superkmer: String::from("AGCAGCTAGCATTTTTGCAGT"),
                minimizer: String::from("AAAAA"),
                is_forward: false,
                start_of_minimizer: 11,
                start_of_superkmer: 0
            }]
        );
    }

    #[test]
    fn test_compute_superkmers2() {
        // same minimizer further away (AAAAA)
        assert_eq!(
            compute_superkmers("AGCAGCTAGCATTTTTGCAGAAAAACC", 16, 5),
            // AGCAGCTAGCATTTTTGCAGAAAA"
            //          CATTTTTGCAGAAAAACC
            vec![
                SuperKmerInfos {
                    superkmer: String::from("AGCAGCTAGCATTTTTGCAGAAAA"),
                    minimizer: String::from("AAAAA"),
                    is_forward: false,
                    start_of_minimizer: 11,
                    start_of_superkmer: 0
                },
                SuperKmerInfos {
                    superkmer: String::from("CATTTTTGCAGAAAAACC"),
                    minimizer: String::from("AAAAA"),
                    is_forward: true,
                    start_of_minimizer: 20,
                    start_of_superkmer: 9
                }
            ]
        );
    }

    #[test]
    fn test_compute_superkmers_sequence_too_short() {
        assert_eq!(compute_superkmers("AGCAGCTAGCATTTT", 16, 5), vec![]);
    }

    #[test]
    fn test_print_() {
        let super_kmer_infos = SuperKmerInfos {
            superkmer: String::from("CATTTTTGCAGAAAAACC"),
            minimizer: String::from("AAAAA"),
            is_forward: true,
            start_of_minimizer: 20,
            start_of_superkmer: 9,
        };
        let s = "SuperKmerInfos { superkmer: \"CATTTTTGCAGAAAAACC\", minimizer: \"AAAAA\", is_forward: true, start_of_minimizer: 20, start_of_superkmer: 9 }";
        assert_eq!(format!("{:?}", super_kmer_infos), s);
    }
}
