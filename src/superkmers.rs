// Get the reverse complement of a DNA sequence
pub fn reverse_complement(seq: &str) -> String {
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
pub fn get_canonical_kmer(kmer: &str) -> (bool, String) {
    let rev_comp = reverse_complement(kmer);
    if kmer < &rev_comp {
        (true, kmer.to_string())
    } else {
        (false, rev_comp.to_string())
    }
}

// OPTIMIZE we don't need so many infos here (the superkmer content itself can be deduced from the original sequence)
#[derive(Debug, PartialEq)]
pub struct SuperKmerInfos {
    pub superkmer: String, // as it would be if the minimizer was canonical in the read
    pub minimizer: String, // in canonical form
    pub was_read_canonical: bool, // was the minimizer in canonical form when reading ?
    pub start_of_minimizer_as_read: usize,
    pub start_of_superkmer_as_read: usize,
}

impl SuperKmerInfos {
    pub fn revcomp_of_superkmer_sequence(&self) -> Self {
        Self {
            superkmer: reverse_complement(&self.superkmer),
            minimizer: self.minimizer.clone(),
            was_read_canonical: self.was_read_canonical,
            start_of_minimizer_as_read: self.start_of_minimizer_as_read,
            start_of_superkmer_as_read: self.start_of_superkmer_as_read,
        }
    }
}

#[derive(PartialEq, Eq)]
struct MinimizerInfos {
    value: String,
    position: usize,
    was_read_canonical: bool,
}

use std::cmp::Ordering;

use crate::brrr_minimizers::MinimizerQueue;
impl PartialOrd for MinimizerInfos {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.value.cmp(&other.value))
    }
}

impl Ord for MinimizerInfos {
    fn cmp(&self, other: &Self) -> Ordering {
        self.value.cmp(&other.value)
    }
}

// Function to compute superkmers
// OPTIMIZE: this is on O(k*sequence.len()). A solution in O(sequence.len()) is possible
pub fn compute_superkmers_linear(sequence: &str, k: usize, m: usize) -> Vec<SuperKmerInfos> {
    let mut superkmers = Vec::with_capacity(sequence.len());
    if sequence.len() < k {
        return superkmers;
    }

    let w = k - m + 1;

    let mut mmers = Vec::with_capacity(sequence.len());
    for position in 0..sequence.len() - m + 1 {
        let (was_read_canonical, value) = get_canonical_kmer(&sequence[position..m + position]);
        mmers.push(MinimizerInfos {
            value,
            position,
            was_read_canonical,
        });
    }

    let mut queue = MinimizerQueue::<_>::new(w);
    let mut minimizers: Vec<(String, usize, bool)> = vec![];

    for (i, input) in mmers.iter().enumerate() {
        if i < w {
            queue.insert(input);
        } else {
            let MinimizerInfos {
                value,
                position,
                was_read_canonical,
            } = queue.get_min();
            minimizers.push((value.clone(), *position, *was_read_canonical));

            queue.insert(input);
        }
    }
    let MinimizerInfos {
        value,
        position,
        was_read_canonical,
    } = queue.get_min();

    minimizers.push((value.clone(), *position, *was_read_canonical));

    let kmer = &sequence[0..k];
    let (current_minimizer, current_relative_position, was_read_canonical) = minimizers[0].clone();
    // superkmer up to here, that will be pushed in the returned vector when the minimizer changes
    let mut superkmerinfos = SuperKmerInfos {
        superkmer: kmer.to_string(), // superkmer will be modified as we see new bases
        minimizer: current_minimizer,
        was_read_canonical,
        start_of_minimizer_as_read: current_relative_position,
        start_of_superkmer_as_read: 0,
    };

    // let mut nb_last_print = 0;
    for i in 1..=sequence.len() - k {
        // if nb_last_print == 0 {
        //     println!("{} (i = {})", i as f32 / sequence.len() as f32 * 100.0, i);
        // }
        // nb_last_print = (nb_last_print + 1) % 10000;

        // for each kmer: if the minimiser the same, append the new base to the superkmer, else create a new superkmer
        let kmer = &sequence[i..i + k];
        let (candidate_minimizer, candidate_absolute_minimizer_position, candidate_is_forward) =
            &minimizers[i];

        // does the minimizer of this kmer have the same value and same position as the minimizer of the previous kmer?
        // let is_same_minimizer = *candidate_minimizer == superkmerinfos.minimizer;
        let is_at_same_position =
            superkmerinfos.start_of_minimizer_as_read == *candidate_absolute_minimizer_position;

        if is_at_same_position {
            // if so, add the new base to the superkmer and keep other data untouched
            superkmerinfos
                .superkmer
                .push(sequence.as_bytes()[i + k - 1] as char);
        } else {
            // otherwise, push the superkmer and create a new one
            if !superkmerinfos.was_read_canonical {
                superkmerinfos = superkmerinfos.revcomp_of_superkmer_sequence();
            }
            superkmers.push(superkmerinfos);
            superkmerinfos = SuperKmerInfos {
                superkmer: kmer.to_string(),
                minimizer: candidate_minimizer.clone(),
                was_read_canonical: *candidate_is_forward,
                start_of_minimizer_as_read: *candidate_absolute_minimizer_position,
                start_of_superkmer_as_read: i,
            };
        }
    }

    if !superkmerinfos.was_read_canonical {
        superkmerinfos = superkmerinfos.revcomp_of_superkmer_sequence();
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
    fn test_compute_superkmers() {
        assert_eq!(
            // ACTGCAAAAATGCTAGCTGCT
            compute_superkmers_linear("AGCAGCTAGCATTTTTGCAGT", 16, 5),
            vec![SuperKmerInfos {
                superkmer: reverse_complement("AGCAGCTAGCATTTTTGCAGT"),
                minimizer: String::from("AAAAA"),
                was_read_canonical: false,
                start_of_minimizer_as_read: 11,
                start_of_superkmer_as_read: 0
            }]
        );
    }
    // encoding and decoding
    #[test]
    fn test_compute_superkmers2() {
        // same minimizer further away (AAAAA)
        assert_eq!(
            compute_superkmers_linear("AGCAGCTAGCATTTTTGCAGAAAAACC", 16, 5),
            // AGCAGCTAGCATTTTTGCAGAAAA"
            //          CATTTTTGCAGAAAAACC
            vec![
                SuperKmerInfos {
                    //TTTTCTGCAAAAATGCTAGCTGCT
                    superkmer: reverse_complement("AGCAGCTAGCATTTTTGCAGAAAA"),
                    minimizer: String::from("AAAAA"),
                    was_read_canonical: false,
                    start_of_minimizer_as_read: 11,
                    start_of_superkmer_as_read: 0
                },
                SuperKmerInfos {
                    superkmer: String::from("CATTTTTGCAGAAAAACC"),
                    minimizer: String::from("AAAAA"),
                    was_read_canonical: true,
                    start_of_minimizer_as_read: 20,
                    start_of_superkmer_as_read: 9
                }
            ]
        );
    }

    // #[test]
    // fn test_compute_superkmers3() {
    //     let superkmers = compute_superkmers_linear("AAGACGCGCCAGCGTCGCATCAGGCGTTGAATGCCGGATGCGCTTCCTGATAAGACGCGCCAGCGTCGCATCAGGCGTTGAATGCCGGATGCGCTT", 40, 31);
    //     for x in &superkmers {
    //         println!("{:?}", x);
    //     }
    //     assert_eq!(superkmers, vec![]);
    // }

    #[test]
    fn test_compute_superkmers_sequence_too_short() {
        assert_eq!(compute_superkmers_linear("AGCAGCTAGCATTTT", 16, 5), vec![]);
    }

    #[test]
    fn test_print_() {
        let super_kmer_infos = SuperKmerInfos {
            superkmer: String::from("CATTTTTGCAGAAAAACC"),
            minimizer: String::from("AAAAA"),
            was_read_canonical: true,
            start_of_minimizer_as_read: 20,
            start_of_superkmer_as_read: 9,
        };
        let s = "SuperKmerInfos { superkmer: \"CATTTTTGCAGAAAAACC\", minimizer: \"AAAAA\", was_read_canonical: true, start_of_minimizer_as_read: 20, start_of_superkmer_as_read: 9 }";
        assert_eq!(format!("{:?}", super_kmer_infos), s);
    }
}
