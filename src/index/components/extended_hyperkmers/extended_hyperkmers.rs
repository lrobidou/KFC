use std::sync::atomic::AtomicUsize;

use crate::{
    buckets::NB_BUCKETS,
    subsequence::{BitPacked, NoBitPacked, Subsequence},
};
use rand::Rng;
use serde::{
    de::{SeqAccess, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};

use std::sync::atomic::Ordering::SeqCst;

use super::array_of_hyperkmer::ArrayOfHyperkmer;

#[derive(PartialEq, Serialize, Deserialize)]
pub struct ParallelExtendedHyperkmers {
    buckets: Vec<ExtendedHyperkmers>,
}

impl ParallelExtendedHyperkmers {
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        let mut buckets = Vec::with_capacity(NB_BUCKETS);
        for _ in 0..NB_BUCKETS {
            buckets.push(ExtendedHyperkmers::new(k, nb_hk_in_a_buffer));
        }
        Self { buckets }
    }

    pub fn get_bucket_from_id(&self, bucket_id: usize) -> &ExtendedHyperkmers {
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

pub struct ExtendedHyperkmers {
    /// the `k` value
    k: usize,
    /// the size in byte taken by a single (extended) hyper kmer
    how_many_u64_for_a_hk: usize,
    /// the number of extended hyperkmer inserted
    nb_inserted: AtomicUsize,
    /// number of extended hyperkmer fitting in a single array
    nb_hk_in_an_array: usize,
    /// the vector of arrays containing extended hyperkmers
    ext_hyperkmers_arrays: boxcar::Vec<ArrayOfHyperkmer>,
    /// the size of each array in u64
    array_size: usize,
}

impl PartialEq for ExtendedHyperkmers {
    fn eq(&self, other: &Self) -> bool {
        let quick_check = self.k == other.k
            && self.how_many_u64_for_a_hk == other.how_many_u64_for_a_hk
            && self.nb_inserted.load(SeqCst) == other.nb_inserted.load(SeqCst)
            && self.nb_hk_in_an_array == other.nb_hk_in_an_array;
        if !quick_check {
            return false;
        }

        // as `self.nb_inserted == other.nb_inserted`, we can do this:
        for index in 0..self.nb_inserted.load(SeqCst) {
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
        state.serialize_field("how_many_u64_for_a_hk", &self.how_many_u64_for_a_hk)?;
        state.serialize_field("nb_inserted", &self.nb_inserted)?;
        state.serialize_field("nb_hk_in_an_array", &self.nb_hk_in_an_array)?;
        state.serialize_field(
            "ext_hyperkmers_arrays",
            &StreamingArrays {
                size: self.how_many_u64_for_a_hk * self.nb_hk_in_an_array,
                arrays: &self.ext_hyperkmers_arrays,
            },
        )?;
        state.end()
    }
}
struct StreamingArrays<'a> {
    size: usize,
    arrays: &'a boxcar::Vec<ArrayOfHyperkmer>,
}

impl<'a> Serialize for StreamingArrays<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeSeq;

        let mut seq = serializer.serialize_seq(Some(self.arrays.count()))?;
        for (_index, array) in self.arrays {
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
                let how_many_u64_for_a_hk = seq
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
                            how_many_u64_for_a_hk * nb_hk_in_an_array,
                        )
                    })
                    .collect();
                let array_size = nb_hk_in_an_array * how_many_u64_for_a_hk;
                Ok(ExtendedHyperkmers {
                    k,
                    how_many_u64_for_a_hk,
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
        let how_many_base_in_a_u64 = 32;
        // TODO simplify
        let how_many_u64_for_a_hk =
            (k - 1) / how_many_base_in_a_u64 + ((k - 1) % how_many_base_in_a_u64 != 0) as usize;
        let array_size = nb_hk_in_an_array * how_many_u64_for_a_hk;

        // pre-allocate one chunk of hyperkmer parts
        let ext_hyperkmers_arrays = boxcar::Vec::new();
        ext_hyperkmers_arrays.push(ArrayOfHyperkmer::new(
            how_many_u64_for_a_hk * nb_hk_in_an_array,
        ));

        Self {
            k,
            ext_hyperkmers_arrays,
            how_many_u64_for_a_hk,
            nb_inserted: AtomicUsize::new(0),
            nb_hk_in_an_array,
            array_size,
        }
    }

    pub fn get_hyperkmer_from_id(&self, id: usize) -> Subsequence<BitPacked> {
        let slice = self.get_slice_from_id(id);
        Subsequence::whole_bitpacked(slice, self.k - 1)
    }

    // SAFETY: no one else should be reading or writing from the memory we are writing to
    unsafe fn dump_hk(&self, id: usize, subseq: Subsequence<NoBitPacked>) {
        debug_assert!(id < self.nb_inserted.load(SeqCst));
        let id_buffer = id / self.nb_hk_in_an_array;
        let pos_in_buffer = id % self.nb_hk_in_an_array; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.how_many_u64_for_a_hk;
        let end = (pos_in_buffer + 1) * (self.how_many_u64_for_a_hk);

        let buffer = &self.ext_hyperkmers_arrays[id_buffer];
        buffer.dump(self.array_size, start, end, subseq);
    }

    fn get_slice_from_id(&self, id: usize) -> &[u64] {
        debug_assert!(id < self.nb_inserted.load(SeqCst));
        let id_buffer = id / self.nb_hk_in_an_array;
        let pos_in_buffer = id % self.nb_hk_in_an_array; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.how_many_u64_for_a_hk;
        let end = (pos_in_buffer + 1) * (self.how_many_u64_for_a_hk);

        let buffer = &self.ext_hyperkmers_arrays[id_buffer];
        &buffer.as_slice(self.array_size)[start..end]
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(&self, new_ext_hyperkmer: &Subsequence<NoBitPacked>) -> usize {
        // let is_full = self.nb_inserted % self.nb_hk_in_an_array == 0;
        let unique_id = self.nb_inserted.fetch_add(1, SeqCst);
        if unique_id % self.nb_hk_in_an_array == 0 {
            // we push to `self.ext_hyperkmers_arrays`
            // while there is some space left (because we pre-allocated it), this signals that we are nearing the end of the capacity
            // TODO annotate that unlikely ???
            self.ext_hyperkmers_arrays.push(ArrayOfHyperkmer::new(
                self.how_many_u64_for_a_hk * self.nb_hk_in_an_array,
            ));
        }

        // dump into memory
        unsafe { self.dump_hk(unique_id, new_ext_hyperkmer.to_canonical()) };

        unique_id
    }

    pub fn get_nb_inserted(&self) -> usize {
        self.nb_inserted.load(SeqCst)
    }
}

impl Drop for ExtendedHyperkmers {
    fn drop(&mut self) {
        for (_index, ext_hyperkmer) in &self.ext_hyperkmers_arrays {
            // SAFETY: no one can access the hyperkmer parts when we are drop it
            unsafe {
                ext_hyperkmer.dealloc(self.k - 1);
            }
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
        let read: Vec<u8> = vec![65, 65, 65];
        let ext_hk = ExtendedHyperkmers::new(k, nb_hk_in_an_array);
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
