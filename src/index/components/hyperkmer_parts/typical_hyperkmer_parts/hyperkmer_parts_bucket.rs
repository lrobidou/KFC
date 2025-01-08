use crate::subsequence::{BitPacked, NoBitPacked, Subsequence};

use super::chunk_of_hyperkmer_parts::ChunkOfHyperkmerParts;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::SeqCst;

use serde::{
    de::{SeqAccess, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};

// TODO rename ?
pub struct HyperkmerPartsBucket {
    /// the `k` value
    k: usize,
    /// the size in byte taken by a single (extended) hyper kmer
    how_many_u64_for_a_hk: usize,
    /// the number of extended hyperkmer inserted
    nb_inserted: AtomicUsize,
    /// number of extended hyperkmer fitting in a single array
    nb_hk_in_an_array: usize,
    /// the vector of arrays containing extended hyperkmers
    ext_hyperkmers_arrays: boxcar::Vec<ChunkOfHyperkmerParts>,
    /// the size of each array in u64
    array_size: usize,
}

impl PartialEq for HyperkmerPartsBucket {
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

impl Serialize for HyperkmerPartsBucket {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("HyperkmerPartsBucket", 5)?;
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
    arrays: &'a boxcar::Vec<ChunkOfHyperkmerParts>,
}

impl Serialize for StreamingArrays<'_> {
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

impl<'de> Deserialize<'de> for HyperkmerPartsBucket {
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

        struct HyperkmerPartsBucketVisitor;

        impl<'de> Visitor<'de> for HyperkmerPartsBucketVisitor {
            type Value = HyperkmerPartsBucket;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("hyperkmer parts bucket")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<HyperkmerPartsBucket, V::Error>
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
                        ChunkOfHyperkmerParts::from_u64_iter(
                            buf,
                            how_many_u64_for_a_hk * nb_hk_in_an_array,
                        )
                    })
                    .collect();
                let array_size = nb_hk_in_an_array * how_many_u64_for_a_hk;
                Ok(HyperkmerPartsBucket {
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
        deserializer.deserialize_struct("HyperkmerPartsBucket", FIELDS, HyperkmerPartsBucketVisitor)
    }
}

impl HyperkmerPartsBucket {
    pub fn new(k: usize, nb_hk_in_an_array: usize) -> Self {
        let how_many_base_in_a_u64 = 32;
        // TODO simplify
        let how_many_u64_for_a_hk =
            (k - 1) / how_many_base_in_a_u64 + ((k - 1) % how_many_base_in_a_u64 != 0) as usize;
        let array_size = nb_hk_in_an_array * how_many_u64_for_a_hk;

        // pre-allocate one chunk of hyperkmer parts
        let ext_hyperkmers_arrays = boxcar::Vec::new();
        ext_hyperkmers_arrays.push(ChunkOfHyperkmerParts::new(
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
            self.ext_hyperkmers_arrays.push(ChunkOfHyperkmerParts::new(
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

impl Drop for HyperkmerPartsBucket {
    fn drop(&mut self) {
        for (_index, ext_hyperkmer) in &self.ext_hyperkmers_arrays {
            // SAFETY: no one can access the hyperkmer parts when we are drop it
            unsafe {
                ext_hyperkmer.dealloc(self.how_many_u64_for_a_hk * self.nb_hk_in_an_array);
            }
        }
    }
}
