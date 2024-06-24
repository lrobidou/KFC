use crate::superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata};

use super::cheap_vec::SimpleVec;

pub struct ExtendedHyperkmers {
    k: usize,
    ext_hyperkmers: Vec<SimpleVec>,
    size_encoded: usize,
}

impl ExtendedHyperkmers {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            ext_hyperkmers: Vec::new(),
            size_encoded: k / 4 + (k % 4 != 0) as usize,
        }
    }

    pub fn get_hyperkmer_from_id(&self, id: usize) -> SubsequenceMetadata<BitPacked> {
        SubsequenceMetadata::whole_bitpacked(
            self.ext_hyperkmers[id].as_slice(self.size_encoded),
            self.k - 1,
        )
    }

    pub fn len(&self) -> usize {
        self.ext_hyperkmers.len()
    }

    /// Adds `new_hyperkmer` in `hyperkmers` and return its index
    /// `new_hyperkmer` does not have to be in canonical form
    pub fn add_new_ext_hyperkmer(
        &mut self,
        new_ext_hyperkmer: &SubsequenceMetadata<NoBitPacked>,
    ) -> usize {
        let id_hyperkmer = self.ext_hyperkmers.len();
        // allocates memory
        let mut new_vec = SimpleVec::new(self.size_encoded);

        // dump into memory
        let dest_slice = new_vec.as_mut_slice(self.size_encoded);
        new_ext_hyperkmer.to_canonical().dump_as_2bits(dest_slice);

        self.ext_hyperkmers.push(new_vec);
        id_hyperkmer
    }
}

impl Drop for ExtendedHyperkmers {
    fn drop(&mut self) {
        for ext_hyperkmer in self.ext_hyperkmers.iter_mut() {
            ext_hyperkmer.dealloc(self.k - 1);
        }
    }
}
