use crate::Index;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

pub fn dump<P: AsRef<Path>>(index: &Index, filename: P) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(filename)?;
    let buffer = BufWriter::with_capacity(9000000, file);
    bincode::serialize_into(buffer, &index)?;
    Ok(())
}

pub fn load<P: AsRef<Path>>(filename: P) -> Result<Index, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let buffer = BufReader::with_capacity(9000000, file);
    let decoded: Index = bincode::deserialize_from(buffer)?;
    Ok(decoded)
}
