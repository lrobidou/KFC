//! This module contains the hyperkmer parts. They are stored in the struct `HyperkmerParts`.
//! There are two kinds of hyperkmer parts: typical ones, and, well, atypical ones.
//!
//! Let's start by describing typical hyperkmer parts. They account for most of the hyperkmer parts. They are always k-1 bases large.
//! They are accessed through `TypicalHyperkmerParts`, which contains `NB_BUCKETS` buckets.
//! Q: Why are the typical parts stored in buckets ?
//! R: So that the total space usage of their storage never doubles when inserting an element.
//! Each bucket contains a vector of chunk of parts. Each chunk contains multiple hyperkmer parts.
//! Q: Why are the parts stored in chunks ? Why not simply a vector of parts ?
//! R: Because allocating a chunk allows doing a single memory allocation for multiple hyperkmer parts.
//!
//! Atypical hyperkmer parts are quite rare: they can appear e.g. 4 times in a genome.
//! They are also larger than k-1.
//! Since they are not the bottleneck, they are simply stored in a vector of tuple (size, content).

mod atypical_hyperkmer_parts;
mod typical_hyperkmer_parts;

use serde::{Deserialize, Serialize};

use atypical_hyperkmer_parts::AtypicalHyperkmerParts;
use typical_hyperkmer_parts::TypicalHyperkmerParts;

use crate::subsequence::{BitPacked, NoBitPacked, Subsequence};

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as likely;
#[cfg(feature = "nightly")]
use core::intrinsics::likely;

use super::HKMetadata;

#[derive(PartialEq, Serialize, Deserialize)]
pub struct HyperkmerParts {
    typical_hyperkmer_parts: TypicalHyperkmerParts,
    atypical_hyperkmer_parts: AtypicalHyperkmerParts,
}

impl HyperkmerParts {
    pub fn new(k: usize, nb_hk_in_a_buffer: usize) -> Self {
        Self {
            typical_hyperkmer_parts: TypicalHyperkmerParts::new(k, nb_hk_in_a_buffer),
            atypical_hyperkmer_parts: AtypicalHyperkmerParts::new(),
        }
    }
    pub fn get_typical_parts(&self) -> &TypicalHyperkmerParts {
        &self.typical_hyperkmer_parts
    }
    pub fn get_atypical_parts(&self) -> &AtypicalHyperkmerParts {
        &self.atypical_hyperkmer_parts
    }
    pub fn add_new_hyperkmer(
        &self,
        is_large: bool,
        seq: &Subsequence<NoBitPacked<'_>>,
    ) -> (usize, usize, bool) {
        if likely(!is_large) {
            let (id_bucket, id_hk) = self.get_typical_parts().add_new_ext_hyperkmer(seq);
            (id_bucket, id_hk, false)
        } else {
            (
                0, // I need an integer here to please the compiler, let's choose 0
                self.get_atypical_parts().add_new_large_hyperkmer(seq),
                true,
            )
        }
    }

    // TODO is this a copy ?
    pub fn get_subsequence_from_metadata<'a>(
        &'a self,
        metadata: &HKMetadata,
    ) -> Subsequence<BitPacked<'a>> {
        let index = metadata.get_index();
        let is_large = metadata.get_is_large();

        if likely(!is_large) {
            self.get_typical_parts()
                .get_bucket_from_id(metadata.get_bucket_id())
                .get_hyperkmer_from_id(index)
        } else {
            let (size, bytes) = self.get_atypical_parts().get(index).unwrap();
            Subsequence::whole_bitpacked(bytes, *size)
        }
    }
}

// #[cfg(any(debug_assertions, test))]
// impl HyperkmerParts {
//     pub fn from_components(
//         typical_hyperkmer_parts: TypicalHyperkmerParts,
//         atypical_hyperkmer_parts: AtypicalHyperkmerParts,
//     ) -> Self {
//         Self {
//             typical_hyperkmer_parts,
//             atypical_hyperkmer_parts,
//         }
//     }
// }

#[cfg(test)]
mod tests {
    use crate::{index::components::HKMetadata, subsequence::Subsequence};

    use super::*;

    #[test]
    fn test_add_large_hyperkmer() {
        let k = 5;
        let nb_hk_in_a_buffer = 13;

        let hyperkmers = HyperkmerParts::new(k, nb_hk_in_a_buffer);

        let large_hyperkmers = hyperkmers.get_atypical_parts();

        let read: Vec<u8> = vec![65, 67, 67, 65, 67, 65, 65, 65, 65, 65, 65];
        let start = 1;
        let end = 8;
        let is_large = true;
        let change_orientation = false;
        let large_hk = Subsequence::new(&read, start, end, !change_orientation);
        let id = large_hyperkmers.add_new_large_hyperkmer(&large_hk);
        assert_eq!(id, 0);
        assert_eq!(large_hyperkmers.len(), 1);
        let metadata = &HKMetadata::new(0, id, start, end, is_large, change_orientation);
        let large_hk_queried = hyperkmers.get_subsequence_from_metadata(metadata);
        assert!(large_hk.equal_bitpacked(&large_hk_queried));
    }
}
