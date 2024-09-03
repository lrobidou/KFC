use std::sync::{Arc, RwLock};

use crate::{
    buckets::{Buckets, NB_BUCKETS},
    subsequence::{BitPacked, NoBitPacked, Subsequence},
};
use rand::Rng;
use serde::{
    de::{SeqAccess, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};

use super::array_of_hyperkmer::ArrayOfHyperkmer;

#[derive(PartialEq, Serialize, Deserialize)]
pub struct ParallelExtendedHyperkmers {
    buckets: Buckets<ExtendedHyperkmers>,
}

impl ParallelExtendedHyperkmers {
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        Self {
            buckets: Buckets::new(|| ExtendedHyperkmers::new(k, nb_hk_in_a_buffer)),
        }
    }

    pub fn get_bucket_from_id_usize(&self, bucket_id: usize) -> Arc<RwLock<ExtendedHyperkmers>> {
        // cloning Arc is cheap
        self.buckets.get_from_id_u64(bucket_id.try_into().unwrap())
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(
        &self,
        new_ext_hyperkmer: &Subsequence<NoBitPacked>,
    ) -> (usize, usize) {
        let mut rng = rand::thread_rng();
        let id_of_chunk: usize = rng.gen_range(0..NB_BUCKETS); // TODO can we use a different number of chunk here

        let extended_hk = self.buckets.get_from_id_usize(id_of_chunk);
        let mut extended_hk = extended_hk.write().unwrap();

        let id_in_this_chunk = extended_hk.add_new_ext_hyperkmer(new_ext_hyperkmer);

        drop(extended_hk);
        (id_of_chunk, id_in_this_chunk)
    }

    pub fn get_nb_inserted(&self) -> usize {
        self.buckets
            .chunks()
            .iter()
            .map(|chunk| chunk.read().unwrap().get_nb_inserted())
            .sum()
    }
}

pub struct ExtendedHyperkmers {
    /// the `k` value
    k: usize,
    /// the size in byte taken by a single (extended) hyper kmer
    byte_size_encoded_hyper_kmer: usize,
    /// the numer of extended hyperkmer inserted
    nb_inserted: usize,
    /// number of extended hyperkmer fitting in a single array
    nb_hk_in_an_array: usize,
    /// the vector of arrays containing extended hyperkmers
    ext_hyperkmers_arrays: Vec<ArrayOfHyperkmer>,
    /// the size of each array
    array_size: usize,
}

impl PartialEq for ExtendedHyperkmers {
    fn eq(&self, other: &Self) -> bool {
        let quick_check = self.k == other.k
            && self.byte_size_encoded_hyper_kmer == other.byte_size_encoded_hyper_kmer
            && self.nb_inserted == other.nb_inserted
            && self.nb_hk_in_an_array == other.nb_hk_in_an_array;
        if !quick_check {
            return false;
        }

        // as `self.nb_inserted == other.nb_inserted`, we can do this:
        for index in 0..self.nb_inserted {
            if self.get_slice_from_id(index) != other.get_slice_from_id(index) {
                return false;
            }
        }
        true
    }
}

impl Serialize for ExtendedHyperkmers {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("ExtendedHyperkmers", 5)?;
        state.serialize_field("k", &self.k)?;
        state.serialize_field(
            "byte_size_encoded_hyper_kmer",
            &self.byte_size_encoded_hyper_kmer,
        )?;
        state.serialize_field("nb_inserted", &self.nb_inserted)?;
        state.serialize_field("nb_hk_in_an_array", &self.nb_hk_in_an_array)?;
        state.serialize_field(
            "ext_hyperkmers_arrays",
            &StreamingArrays {
                size: self.byte_size_encoded_hyper_kmer * self.nb_hk_in_an_array,
                arrays: &self.ext_hyperkmers_arrays,
            },
        )?;
        state.end()
    }
}
struct StreamingArrays<'a> {
    size: usize,
    arrays: &'a [ArrayOfHyperkmer],
}

impl<'a> Serialize for StreamingArrays<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeSeq;

        let mut seq = serializer.serialize_seq(Some(self.arrays.len()))?;
        for array in self.arrays {
            seq.serialize_element(array.as_u64_slice(self.size))?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for ExtendedHyperkmers {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "snake_case")]
        enum Field {
            K,
            SizeEncoded,
            NbInserted,
            NbHkInABuffer,
            ExtHyperkmersBuffers,
        }

        struct ExtendedHyperkmersVisitor;

        impl<'de> Visitor<'de> for ExtendedHyperkmersVisitor {
            type Value = ExtendedHyperkmers;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("extended hyperkmers")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<ExtendedHyperkmers, V::Error>
            where
                V: SeqAccess<'de>,
            {
                let k = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
                let byte_size_encoded_hyper_kmer = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(1, &self))?;
                let nb_inserted = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(2, &self))?;
                let nb_hk_in_an_array = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(3, &self))?;
                let arrays: Vec<Vec<u64>> = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(4, &self))?;
                let ext_hyperkmers_arrays = arrays
                    .into_iter()
                    .map(|buf| {
                        ArrayOfHyperkmer::from_u64_iter(
                            buf,
                            byte_size_encoded_hyper_kmer * nb_hk_in_an_array,
                        )
                    })
                    .collect();
                let array_size = (nb_hk_in_an_array * byte_size_encoded_hyper_kmer) / 8
                    + ((nb_hk_in_an_array * byte_size_encoded_hyper_kmer) % 8 != 0) as usize;
                Ok(ExtendedHyperkmers {
                    k,
                    byte_size_encoded_hyper_kmer,
                    nb_inserted,
                    nb_hk_in_an_array,
                    ext_hyperkmers_arrays,
                    array_size,
                })
            }
        }

        const FIELDS: &[&str] = &[
            "k",
            "byte_size_encoded_hyper_kmer",
            "nb_inserted",
            "nb_hk_in_a_buffer",
            "ext_hyperkmers_buffers",
        ];
        deserializer.deserialize_struct("ExtendedHyperkmers", FIELDS, ExtendedHyperkmersVisitor)
    }
}

impl ExtendedHyperkmers {
    pub fn new(k: usize, nb_hk_in_an_array: usize) -> Self {
        let byte_size_encoded_hyper_kmer = (k - 1) / 4 + ((k - 1) % 4 != 0) as usize;
        let array_size = (nb_hk_in_an_array * byte_size_encoded_hyper_kmer) / 8
            + ((nb_hk_in_an_array * byte_size_encoded_hyper_kmer) % 8 != 0) as usize;
        Self {
            k,
            ext_hyperkmers_arrays: Vec::new(),
            byte_size_encoded_hyper_kmer,
            nb_inserted: 0,
            nb_hk_in_an_array,
            array_size,
        }
    }

    pub fn get_hyperkmer_from_id(&self, id: usize) -> Subsequence<BitPacked> {
        let slice = self.get_slice_from_id(id);
        Subsequence::whole_bitpacked(slice, self.k - 1)
    }

    pub fn len(&self) -> usize {
        self.nb_inserted
    }

    fn dump_hk(&mut self, id: usize, subseq: Subsequence<NoBitPacked>) {
        debug_assert!(id < self.nb_inserted);
        let id_buffer = id / self.nb_hk_in_an_array;
        let pos_in_buffer = id % self.nb_hk_in_an_array; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.byte_size_encoded_hyper_kmer;
        let end = (pos_in_buffer + 1) * (self.byte_size_encoded_hyper_kmer);

        let buffer = &mut self.ext_hyperkmers_arrays[id_buffer];
        buffer.dump(self.array_size, start, end, subseq);
    }

    fn get_slice_from_id(&self, id: usize) -> &[u8] {
        debug_assert!(id < self.nb_inserted);
        let id_buffer = id / self.nb_hk_in_an_array;
        let pos_in_buffer = id % self.nb_hk_in_an_array; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.byte_size_encoded_hyper_kmer;
        let end = (pos_in_buffer + 1) * (self.byte_size_encoded_hyper_kmer);

        let buffer = &self.ext_hyperkmers_arrays[id_buffer];
        &buffer.as_slice(self.array_size)[start..end]
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(&mut self, new_ext_hyperkmer: &Subsequence<NoBitPacked>) -> usize {
        let is_full = self.nb_inserted % self.nb_hk_in_an_array == 0;

        // allocates memory
        if is_full {
            self.ext_hyperkmers_arrays.push(ArrayOfHyperkmer::new(
                self.byte_size_encoded_hyper_kmer * self.nb_hk_in_an_array,
            ));
        }
        let id_hyperkmer = self.len();
        self.nb_inserted += 1;

        // dump into memory
        self.dump_hk(id_hyperkmer, new_ext_hyperkmer.to_canonical());

        id_hyperkmer
    }

    pub fn get_nb_inserted(&self) -> usize {
        self.nb_inserted
    }
}

impl Drop for ExtendedHyperkmers {
    fn drop(&mut self) {
        for ext_hyperkmer in self.ext_hyperkmers_arrays.iter_mut() {
            ext_hyperkmer.dealloc(self.k - 1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extended_hk_get_nb_inserted() {
        let k = 31;
        let nb_hk_in_an_array = 7;
        let read: Vec<u8> = vec![3, 8, 7];
        let mut ext_hk = ExtendedHyperkmers::new(k, nb_hk_in_an_array);
        for _ in 0..10 {
            ext_hk.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        assert_eq!(ext_hk.get_nb_inserted(), 10);
    }

    #[test]
    fn test_parallel_extended_hk_get_nb_inserted() {
        let k = 31;
        let nb_hk_in_a_buffer = 7;

        let ext_hks = ParallelExtendedHyperkmers::new(k, nb_hk_in_a_buffer);

        let ext_hk_6 = ext_hks.get_bucket_from_id_usize(6);
        let mut ext_hk_6 = ext_hk_6.write().unwrap();
        let ext_hk_8 = ext_hks.get_bucket_from_id_usize(8);
        let mut ext_hk_8 = ext_hk_8.write().unwrap();
        let ext_hk_9: Arc<RwLock<ExtendedHyperkmers>> = ext_hks.get_bucket_from_id_usize(9);
        let mut ext_hk_9 = ext_hk_9.write().unwrap();
        let read: Vec<u8> = vec![3, 8, 7];
        for _ in 0..10 {
            ext_hk_6.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        let read: Vec<u8> = vec![3, 8, 7, 9, 0, 7];
        for _ in 0..10 {
            ext_hk_8.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        let read: Vec<u8> = vec![0, 3, 8, 7, 0];
        for _ in 0..10 {
            ext_hk_9.add_new_ext_hyperkmer(&Subsequence::new(&read, 0, read.len(), true));
        }
        drop(ext_hk_6);
        drop(ext_hk_8);
        drop(ext_hk_9);
        assert_eq!(ext_hks.get_nb_inserted(), 30);
    }
}
