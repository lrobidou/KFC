use crate::{Count, HashSuperKmer, Minimizer, SuperKmerInfos};
use mashmap::MashMap;
use xxhash_rust::const_xxh3::xxh3_64;

pub struct SuperKmerCounts {
    data: MashMap<Minimizer, (HashSuperKmer, Count)>,
}

impl SuperKmerCounts {
    pub fn new() -> Self {
        Self {
            data: MashMap::new(),
        }
    }
    /// Retrieve the count of `superkmer` in `sk_count`.
    /// the superkmer and its minimizer should be canonical.
    /// Returns 0 if the superkmer is not found.
    pub fn get_count_superkmer(&self, superkmer: &SuperKmerInfos) -> Count {
        let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
        for (hash, count) in self.data.get_mut_iter(&superkmer.minimizer) {
            if *hash == superkmer_hash {
                return *count;
            }
        }
        // superkmer not found, count is 0
        0
    }

    /// Add 1 to the count of `superkmer` in `sk_count`.
    /// If `superkmer` is not in `sk_count`, add it.
    /// Returns the new count.
    pub fn increase_count_superkmer(&mut self, superkmer: &SuperKmerInfos) -> Count {
        let superkmer_hash = xxh3_64(superkmer.superkmer.as_bytes());
        for (super_kmer_hash, super_kmer_count) in self.data.get_mut_iter(&superkmer.minimizer) {
            if *super_kmer_hash == superkmer_hash {
                let new_count = super_kmer_count.saturating_add(1);
                *super_kmer_count = new_count;
                return new_count;
            }
        }
        // superkmer not found, insert it, count is 1
        &self
            .data
            .insert(superkmer.minimizer.clone(), (superkmer_hash, 1));
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_count_superkmer() {
        let sk0 = SuperKmerInfos {
            superkmer: "ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGG".into(),
            minimizer: "AATCGTCATCCGAAACGTCA".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 7,
            start_of_superkmer_as_read: 0,
        };
        let sk1 = SuperKmerInfos {
            superkmer: "GACGTTTCGGATGACGATTGTACGTGACGGTG".into(),
            minimizer: "ACAATCGTCATCCGAAACGT".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 9,
            start_of_superkmer_as_read: 8,
        };
        let sk2 = SuperKmerInfos {
            superkmer: "CGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATG".into(),
            minimizer: "ACCGTCACGTACAATCGTCA".into(),
            was_read_canonical: false,
            start_of_minimizer_as_read: 19,
            start_of_superkmer_as_read: 10,
        };
        let sk3 = SuperKmerInfos {
            superkmer: "GACGATTGTACGTGACGGTGCGTCCGGATGAC".into(),
            minimizer: "ACGATTGTACGTGACGGTGC".into(),
            was_read_canonical: true,
            start_of_minimizer_as_read: 21,
            start_of_superkmer_as_read: 20,
        };

        let mut sk_count = SuperKmerCounts::new();

        for sk in vec![sk0, sk1, sk2, sk3] {
            assert_eq!(sk_count.get_count_superkmer(&sk), 0);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 1);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 2);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 3);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 4);
            assert_eq!(sk_count.get_count_superkmer(&sk), 4);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 5);
            assert_eq!(sk_count.increase_count_superkmer(&sk), 6);
            assert_eq!(sk_count.get_count_superkmer(&sk), 6);
        }
    }
}
