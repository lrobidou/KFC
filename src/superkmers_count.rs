use crate::{Count, HashSuperKmer, Minimizer, Superkmer};
use mashmap::MashMap;

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
    pub fn get_count_superkmer(&self, superkmer: &Superkmer) -> Count {
        let superkmer_hash = superkmer.hash_superkmer();
        for (hash, count) in self.data.get_mut_iter(&superkmer.get_minimizer()) {
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
    pub fn increase_count_superkmer(&mut self, superkmer: &Superkmer) -> Count {
        let superkmer_hash = superkmer.hash_superkmer();
        for (super_kmer_hash, super_kmer_count) in
            self.data.get_mut_iter(&superkmer.get_minimizer())
        {
            if *super_kmer_hash == superkmer_hash {
                let new_count = super_kmer_count.saturating_add(1);
                *super_kmer_count = new_count;
                return new_count;
            }
        }
        // superkmer not found, insert it, count is 1
        self.data
            .insert(superkmer.get_minimizer(), (superkmer_hash, 1));
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_count_superkmer() {
        let m = 10;
        // read: ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATGAC
        //       ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGG                 (AATCGTCATCCGAAACGTCA, 7)
        //               GACGTTTCGGATGACGATTGTACGTGACGGTG               (ACAATCGTCATCCGAAACGT, 9)
        //                 CGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATG     (ACCGTCACGTACAATCGTCA, 19)
        //                           GACGATTGTACGTGACGGTGCGTCCGGATGAC   (ACGATTGTACGTGACGGTGC, 21)
        let read = "ACGTACGTGACGTTTCGGATGACGATTGTACGTGACGGTGCGTCCGGATGAC";
        let sk0 = Superkmer::new(read, 7, 7 + m, 0, 38, false);
        let sk1 = Superkmer::new(read, 9, 9 + m, 9, 9 + 32, false);
        let sk2 = Superkmer::new(read, 19, 19 + m, 10, 10 + 40, false);
        let sk3 = Superkmer::new(read, 21, 21 + m, 20, 20 + 32, true);

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
