use std::{fs::File, path::Path};

use crate::Index;

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn plain_text<P: AsRef<Path>>(index: &Index, output_file: P) {
    // use std::io::Write;
    // let mut file = File::create(output_file).unwrap();

    todo!();
    // // Iterate over the vector and write each line to the file
    // for i in 0..ext_hyperkmers.len() {
    //     let ext_hk = ext_hyperkmers.get_hyperkmer_from_id(i);
    //     let line = ext_hk.to_string();
    //     writeln!(file, "{line}").unwrap();
    // }
}
