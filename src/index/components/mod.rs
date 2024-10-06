mod extended_hyperkmers;
mod hyperkmers_counts;
mod superkmers_count;

pub use extended_hyperkmers::{ExtendedHyperkmers, ParallelExtendedHyperkmers};
pub use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
pub use superkmers_count::SuperKmerCounts;
pub type LargeExtendedHyperkmer = (usize, Vec<u64>);

use crate::subsequence::{BitPacked, NoBitPacked, Subsequence};

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as likely;
#[cfg(feature = "nightly")]
use core::intrinsics::likely;
use std::sync::{Arc, RwLock};

pub fn get_subsequence_from_metadata<'a>(
    hyperkmers: &'a ExtendedHyperkmers,
    large_hyperkmers: &'a [LargeExtendedHyperkmer],
    metadata: &HKMetadata,
) -> Subsequence<BitPacked<'a>> {
    let index = metadata.get_index();
    let is_large = metadata.get_is_large();
    if likely(!is_large) {
        hyperkmers.get_hyperkmer_from_id(index)
    } else {
        let (size, bytes) = &large_hyperkmers[index];
        Subsequence::whole_bitpacked(bytes, *size)
    }
}

// Write locks on a hyperkmer
pub fn add_new_hyperkmer(
    is_large: bool,
    seq: &Subsequence<NoBitPacked<'_>>,
    hyperkmers: &ParallelExtendedHyperkmers,
    large_hyperkmers: &Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
) -> (usize, usize, bool) {
    if likely(!is_large) {
        let (id_bucket, id_hk) = hyperkmers.add_new_ext_hyperkmer(seq);
        (id_bucket, id_hk, false)
    } else {
        let mut large_hyperkmers = large_hyperkmers.write().unwrap();
        (
            0, // I need an integer here to please the compiler, let's choose 0
            add_new_large_hyperkmer(&mut large_hyperkmers, seq),
            true,
        )
    }
}

/// Adds `new_ext_hyperkmer` into `large_hyperkmers`.
///
/// Returns the position of `new_ext_hyperkmer` in `large_hyperkmers`
pub fn add_new_large_hyperkmer(
    large_hyperkmers: &mut Vec<LargeExtendedHyperkmer>,
    new_ext_hyperkmer: &Subsequence<NoBitPacked>,
) -> usize {
    let nb_base = new_ext_hyperkmer.len();

    let size_needed_for_slice = nb_base / 32 + (nb_base % 32 != 0) as usize;

    let mut dest_slice = vec![0; size_needed_for_slice];
    new_ext_hyperkmer
        .to_canonical()
        .dump_as_2bits(&mut dest_slice);

    large_hyperkmers.push((nb_base, dest_slice));
    large_hyperkmers.len() - 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_large_hyperkmer() {
        let k = 5;
        let nb_hk_in_a_buffer = 13;
        let hyperkmers = ExtendedHyperkmers::new(k, nb_hk_in_a_buffer);
        let mut large_hyperkmers = Vec::<LargeExtendedHyperkmer>::new();

        let read: Vec<u8> = vec![65, 67, 67, 65, 67, 65, 65, 65, 65, 65, 65];
        let start = 1;
        let end = 8;
        let is_large = true;
        let change_orientation = false;
        let large_hk = Subsequence::new(&read, start, end, !change_orientation);
        let id = add_new_large_hyperkmer(&mut large_hyperkmers, &large_hk);
        assert_eq!(id, 0);
        assert_eq!(large_hyperkmers.len(), 1);
        let metadata = &HKMetadata::new(0, id, start, end, is_large, change_orientation);
        let large_hk_queried =
            get_subsequence_from_metadata(&hyperkmers, &large_hyperkmers, metadata);
        assert!(large_hk.equal_bitpacked(&large_hk_queried));
    }
}
