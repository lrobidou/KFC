use std::{fs::File, path::Path};

use crate::extended_hyperkmers::ExtendedHyperkmers;

pub fn plain_text<P: AsRef<Path>>(ext_hyperkmers: &ExtendedHyperkmers, output_file: P) {
    use std::io::Write;
    let mut file = File::create(output_file).unwrap();

    // Iterate over the vector and write each line to the file
    for i in 0..ext_hyperkmers.len() {
        let ext_hk = ext_hyperkmers.get_hyperkmer_from_id(i);
        let line = ext_hk.to_string();
        writeln!(file, "{line}").unwrap();
    }
}
