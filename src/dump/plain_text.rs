use std::{fs::File, io::BufWriter, path::Path};

use crate::Index;

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn plain_text<P: AsRef<Path>>(index: &Index, output_file: P) {
    use std::io::Write;
    let file = File::create(output_file).unwrap();
    let mut buffer = BufWriter::new(file);
    for (kmer, count) in index.iter_kmers().flatten() {
        writeln!(buffer, "{}\t{}", kmer, count).unwrap();
    }
}
