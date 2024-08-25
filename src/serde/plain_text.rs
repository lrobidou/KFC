use std::{fs::File, io::BufWriter, path::Path};

use serde::Serialize;

use crate::{index::FullIndexTrait, Index};

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn plain_text<P: AsRef<Path>, FI: FullIndexTrait + Serialize + Send + Sync>(
    index: &Index<FI>,
    output_file: P,
) {
    // use std::io::Write;
    // let file = File::create(output_file).unwrap();
    // let mut buffer = BufWriter::with_capacity(100_000, file);
    // for (kmer, count) in index.iter_kmers().flatten() {
    //     buffer
    //         .write_all(&kmer)
    //         .expect("could not write to the file");
    //     writeln!(buffer, "\t{}", count).unwrap();
    // }
    index.par_write_kmers(output_file);
}
