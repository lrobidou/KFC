use crate::superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata};

use super::cheap_vec::SimpleVec;

pub struct ExtendedHyperkmers {
    k: usize,
    ext_hyperkmers_buffers: Vec<SimpleVec>,
    size_encoded: usize,
    nb_inserted: usize,
    nb_hk_in_a_buffer: usize,
}

impl ExtendedHyperkmers {
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        Self {
            k,
            ext_hyperkmers_buffers: Vec::new(),
            size_encoded: (k - 1) / 4 + ((k - 1) % 4 != 0) as usize,
            nb_inserted: 0,
            nb_hk_in_a_buffer,
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
        let pos_in_buffer = id % self.nb_hk_in_a_buffer;

        let buffer = &mut self.ext_hyperkmers_buffers[id_buffer];
        &mut buffer.as_mut_slice(self.size_encoded * self.nb_hk_in_a_buffer)
            [pos_in_buffer * self.size_encoded..(pos_in_buffer + 1) * self.size_encoded]
    }

    fn get_slice_from_id(&self, id: usize) -> &[u8] {
        debug_assert!(id < self.nb_inserted);
        let id_buffer = id / self.nb_hk_in_a_buffer;
        let pos_in_buffer = id % self.nb_hk_in_a_buffer;

        let buffer = &self.ext_hyperkmers_buffers[id_buffer];
        &buffer.as_slice(self.size_encoded * self.nb_hk_in_a_buffer)
            [pos_in_buffer * self.size_encoded..(pos_in_buffer + 1) * (self.size_encoded)]
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
            self.ext_hyperkmers_buffers
                .push(SimpleVec::new(self.size_encoded * self.nb_hk_in_a_buffer));
        }
        let id_hyperkmer = self.len();
        self.nb_inserted += 1;

        // dump into memory
        let dest_slice = self.get_mut_slice_from_id(id_hyperkmer);
        new_ext_hyperkmer.to_canonical().dump_as_2bits(dest_slice);

        id_hyperkmer
    }
}

impl Drop for ExtendedHyperkmers {
    fn drop(&mut self) {
        for ext_hyperkmer in self.ext_hyperkmers_buffers.iter_mut() {
            ext_hyperkmer.dealloc(self.k - 1);
        }
    }
}
