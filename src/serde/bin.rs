use serde::{Deserialize, Serialize};

use crate::index::FullIndexTrait;
use crate::Index;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

pub fn dump<P: AsRef<Path>, FI: FullIndexTrait + Serialize + Send + Sync + Serialize>(
    index: &Index<FI>,
    filename: P,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(filename)?;
    let buffer = BufWriter::with_capacity(9000000, file);
    bincode::serialize_into(buffer, &index)?;
    Ok(())
}

pub fn load<
    P: AsRef<Path>,
    FI: FullIndexTrait + Send + Sync + Serialize + for<'a> Deserialize<'a>,
>(
    filename: P,
) -> Result<Index<FI>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let buffer = BufReader::with_capacity(9000000, file);
    let decoded: Index<FI> = bincode::deserialize_from(buffer)?;
    Ok(decoded)
}
