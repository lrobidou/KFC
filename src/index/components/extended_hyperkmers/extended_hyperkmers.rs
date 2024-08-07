use serde::{
    de::{SeqAccess, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};

use crate::superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata};

use super::cheap_vec::SimpleVec;
use crate::macros::debug_print as p;

pub struct ExtendedHyperkmers {
    /// the `k` value
    k: usize,
    /// the size in byte taken by a single (extended) hyper kmer
    byte_size_encoded_hyper_kmer: usize,
    /// the numer of extended hyperkmer inserted
    nb_inserted: usize,
    /// number of extended hyperkmer fitting in a single buffer
    nb_hk_in_a_buffer: usize,
    /// the vector of buffers containing extended hyperkmers
    ext_hyperkmers_buffers: Vec<SimpleVec>,
    /// the size of each buffer
    buffer_size: usize,
}

impl PartialEq for ExtendedHyperkmers {
    fn eq(&self, other: &Self) -> bool {
        let quick_check = self.k == other.k
            && self.byte_size_encoded_hyper_kmer == other.byte_size_encoded_hyper_kmer
            && self.nb_inserted == other.nb_inserted
            && self.nb_hk_in_a_buffer == other.nb_hk_in_a_buffer;
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
        let mut state = serializer.serialize_struct("ExtendedHyperkmers", 0)?;
        state.serialize_field("k", &self.k)?;
        state.serialize_field(
            "byte_size_encoded_hyper_kmer",
            &self.byte_size_encoded_hyper_kmer,
        )?;
        state.serialize_field("nb_inserted", &self.nb_inserted)?;
        state.serialize_field("nb_hk_in_a_buffer", &self.nb_hk_in_a_buffer)?;
        state.serialize_field(
            "ext_hyperkmers_buffers",
            &StreamingBuffers {
                size: self.byte_size_encoded_hyper_kmer * self.nb_hk_in_a_buffer,
                buffers: &self.ext_hyperkmers_buffers,
            },
        )?;
        state.end()
    }
}
struct StreamingBuffers<'a> {
    size: usize,
    buffers: &'a [SimpleVec],
}

impl<'a> Serialize for StreamingBuffers<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeSeq;

        let mut seq = serializer.serialize_seq(Some(self.buffers.len()))?;
        for buffer in self.buffers {
            seq.serialize_element(buffer.as_u64_slice(self.size))?;
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
                let nb_hk_in_a_buffer = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(3, &self))?;
                let buffers: Vec<Vec<u64>> = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(4, &self))?;
                let ext_hyperkmers_buffers = buffers
                    .into_iter()
                    .map(|buf| {
                        SimpleVec::from_u64_iter(
                            buf,
                            byte_size_encoded_hyper_kmer * nb_hk_in_a_buffer,
                        )
                    })
                    .collect();
                let buffer_size = (nb_hk_in_a_buffer * byte_size_encoded_hyper_kmer) / 8
                    + ((nb_hk_in_a_buffer * byte_size_encoded_hyper_kmer) % 8 != 0) as usize;
                Ok(ExtendedHyperkmers {
                    k,
                    byte_size_encoded_hyper_kmer,
                    nb_inserted,
                    nb_hk_in_a_buffer,
                    ext_hyperkmers_buffers,
                    buffer_size,
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
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        let byte_size_encoded_hyper_kmer = (k - 1) / 4 + ((k - 1) % 4 != 0) as usize;
        let buffer_size = (nb_hk_in_a_buffer * byte_size_encoded_hyper_kmer) / 8
            + ((nb_hk_in_a_buffer * byte_size_encoded_hyper_kmer) % 8 != 0) as usize;
        Self {
            k,
            ext_hyperkmers_buffers: Vec::new(),
            byte_size_encoded_hyper_kmer,
            nb_inserted: 0,
            nb_hk_in_a_buffer,
            buffer_size,
        }
    }

    pub fn get_hyperkmer_from_id(&self, id: usize) -> SubsequenceMetadata<BitPacked> {
        let slice = self.get_slice_from_id(id);
        SubsequenceMetadata::whole_bitpacked(slice, self.k - 1)
    }

    pub fn len(&self) -> usize {
        self.nb_inserted
    }

    fn get_mut_slice_from_id(&mut self, id: usize) -> &mut [u8] {
        debug_assert!(id < self.nb_inserted);
        let id_buffer = id / self.nb_hk_in_a_buffer;
        let pos_in_buffer = id % self.nb_hk_in_a_buffer; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.byte_size_encoded_hyper_kmer;
        let end = (pos_in_buffer + 1) * (self.byte_size_encoded_hyper_kmer);

        let buffer = &mut self.ext_hyperkmers_buffers[id_buffer];
        &mut buffer.as_mut_slice(self.buffer_size)[start..end]
    }

    fn get_slice_from_id(&self, id: usize) -> &[u8] {
        debug_assert!(id < self.nb_inserted);
        let id_buffer = id / self.nb_hk_in_a_buffer;
        let pos_in_buffer = id % self.nb_hk_in_a_buffer; // which hyperkmer is it from the buffer `id_buffer`?
        let start = pos_in_buffer * self.byte_size_encoded_hyper_kmer;
        let end = (pos_in_buffer + 1) * (self.byte_size_encoded_hyper_kmer);

        let buffer = &self.ext_hyperkmers_buffers[id_buffer];
        &buffer.as_slice(self.buffer_size)[start..end]
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(
        &mut self,
        new_ext_hyperkmer: &SubsequenceMetadata<NoBitPacked>,
    ) -> usize {
        let is_full = self.nb_inserted % self.nb_hk_in_a_buffer == 0;

        // allocates memory
        if is_full {
            self.ext_hyperkmers_buffers.push(SimpleVec::new(
                self.byte_size_encoded_hyper_kmer * self.nb_hk_in_a_buffer,
            ));
        }
        let id_hyperkmer = self.len();
        self.nb_inserted += 1;

        // dump into memory
        let dest_slice = self.get_mut_slice_from_id(id_hyperkmer);
        new_ext_hyperkmer.to_canonical().dump_as_2bits(dest_slice);

        id_hyperkmer
    }

    pub fn get_nb_inserted(&self) -> usize {
        self.nb_inserted
    }
}

impl Drop for ExtendedHyperkmers {
    fn drop(&mut self) {
        for ext_hyperkmer in self.ext_hyperkmers_buffers.iter_mut() {
            ext_hyperkmer.dealloc(self.k - 1);
        }
    }
}