use std::path::Path;

use serde::Serialize;

use crate::{index::FullIndexTrait, Count, Index};

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn plain_text<P: AsRef<Path>, FI: FullIndexTrait + Serialize + Send + Sync>(
    index: &Index<FI>,
    kmer_threshold: Count,
    output_file: P,
) {
    index.par_write_kmers(kmer_threshold, output_file);
}
