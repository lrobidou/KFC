mod extended_hyperkmers;
mod hyperkmers_counts;
mod superkmers_count;

pub use extended_hyperkmers::ExtendedHyperkmers;
pub use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
pub use superkmers_count::SuperKmerCounts;

use crate::superkmer::{BitPacked, NoBitPacked, SubsequenceMetadata};

pub fn get_subsequence_from_metadata<'a>(
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a Vec<(usize, Vec<u8>)>,
    metadata: &HKMetadata,
) -> SubsequenceMetadata<'a, BitPacked> {
    let index = metadata.get_index();
    let is_large = metadata.get_is_large();
    // TODO likely
    if !is_large {
        hyperkmers.get_hyperkmer_from_id(index)
    } else {
        let (size, bytes) = &large_hyperkmers[index];
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
