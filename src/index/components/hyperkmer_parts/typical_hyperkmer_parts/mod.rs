mod chunk_of_hyperkmer_parts;
mod hyperkmer_parts_bucket;

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{
    buckets::NB_BUCKETS,
    subsequence::{NoBitPacked, Subsequence},
};

pub use hyperkmer_parts_bucket::HyperkmerPartsBucket;

#[derive(PartialEq, Serialize, Deserialize)]
pub struct TypicalHyperkmerParts {
    buckets: Vec<HyperkmerPartsBucket>,
}

impl TypicalHyperkmerParts {
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        let mut buckets = Vec::with_capacity(NB_BUCKETS);
        for _ in 0..NB_BUCKETS {
            buckets.push(HyperkmerPartsBucket::new(k, nb_hk_in_a_buffer));
        }
        Self { buckets }
    }

    pub fn get_bucket_from_id(&self, bucket_id: usize) -> &HyperkmerPartsBucket {
        &self.buckets[bucket_id]
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(
        &self,
        new_ext_hyperkmer: &Subsequence<NoBitPacked>,
    ) -> (usize, usize) {
        let mut rng = rand::thread_rng();
        let id_of_chunk: usize = rng.gen_range(0..NB_BUCKETS); // TODO can we use a different number of chunk here

        let extended_hk = self.get_bucket_from_id(id_of_chunk);

        let id_in_this_chunk = extended_hk.add_new_ext_hyperkmer(new_ext_hyperkmer);

        (id_of_chunk, id_in_this_chunk)
    }

    pub fn get_nb_inserted(&self) -> usize {
        self.buckets
            .iter()
            .map(|chunk| chunk.get_nb_inserted())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extended_hk_get_nb_inserted() {
        let k = 31;
        let nb_hk_in_an_array = 7;
        let read: Vec<u8> = vec![65, 65, 65];
        let ext_hk = HyperkmerPartsBucket::new(k, nb_hk_in_an_array);
        for _ in 0..10 {
            ext_hk.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        assert_eq!(ext_hk.get_nb_inserted(), 10);
    }

    #[test]
    fn test_parallel_extended_hk_get_nb_inserted() {
        let k = 31;
        let nb_hk_in_a_buffer = 7;

        let ext_hks = TypicalHyperkmerParts::new(k, nb_hk_in_a_buffer);

        // caution: to prevent deadlock, the locks are acquired in ascending order.
        // so this test might deadlock if the number of buckets is smaller than 9
        assert!(ext_hks.buckets.len() > 9);

        let ext_hk_6 = ext_hks.get_bucket_from_id(6);
        let ext_hk_8 = ext_hks.get_bucket_from_id(8);
        let ext_hk_9 = ext_hks.get_bucket_from_id(9);
        let read: Vec<u8> = vec![b'N', b'A', b'A'];
        for _ in 0..10 {
            ext_hk_6.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        let read: Vec<u8> = vec![b'N', b'A', b'A', b'N', b'A', b'A'];
        for _ in 0..10 {
            ext_hk_8.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        let read: Vec<u8> = vec![b'A', b'A', b'N', b'A'];
        for _ in 0..10 {
            ext_hk_9.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        assert_eq!(ext_hks.get_nb_inserted(), 30);
    }
}
