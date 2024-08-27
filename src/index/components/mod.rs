mod extended_hyperkmers;
mod hyperkmers_counts;
mod superkmers_count;

pub use extended_hyperkmers::ParallelExtendedHyperkmers;
pub use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
pub use superkmers_count::SuperKmerCounts;

use crate::superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata};

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as likely;
#[cfg(feature = "nightly")]
use core::intrinsics::likely;

pub fn get_subsequence_from_metadata<'a>(
    hyperkmers: &'a ParallelExtendedHyperkmers,
    large_hyperkmers: &'a [(usize, Vec<u8>)],
    metadata: &HKMetadata,
) -> SubsequenceMetadata<BitPacked> {
    let index = metadata.get_index();
    let is_large = metadata.get_is_large();
    if likely(!is_large) {
        hyperkmers.get_hyperkmer_from_id(index)
    } else {
        let (size, bytes) = &large_hyperkmers[index];
        let bytes = bytes.clone();
        let subseq = crate::superkmer::SubsequenceMetadata::whole_bitpacked(bytes, *size);
        subseq
    }
}

pub fn add_new_large_hyperkmer(
    large_hyperkmers: &mut Vec<(usize, Vec<u8>)>,
    new_ext_hyperkmer: &SubsequenceMetadata<NoBitPacked>,
) -> usize {
    let nb_base = new_ext_hyperkmer.len();

    let size_needed_for_slice = nb_base / 4 + (nb_base % 4 != 0) as usize;

    let mut dest_slice = vec![0; size_needed_for_slice];
    new_ext_hyperkmer
        .to_canonical()
        .dump_as_2bits(&mut dest_slice);

    large_hyperkmers.push((nb_base, dest_slice));
    large_hyperkmers.len() - 1
}
