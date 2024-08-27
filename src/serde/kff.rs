use crate::index::FullIndexTrait;
use crate::{Count, Index};
use bitvec::view::BitView;
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use kff::Kff;
use log::error;
use serde::Serialize;
use std::{fs::File, io::BufWriter, path::Path};
use thiserror::Error;

pub const ENCODING: u8 = 0b00011011;

/// Error related to manipulation of KFF files
#[derive(Error, Debug)]
pub enum KFFError {
    #[error("Failed to cast values to a KFF-compatible format")]
    FailedCast,
    #[error("Failed to write values on a KFF file")]
    WriteValues,
    #[error("Missing field")]
    MissingField,
    #[error("Could not write the minimizer section")]
    WriteMinimizer,
    #[error("Could not finalize the kff file")]
    FinalizeFailure,
    #[error("Invalid header")]
    InvalidHeader,
    #[error("Could not create the kff file")]
    CreationFailure,
}

/// macro to write a function that logs an error and returns `KFFError::FailedCast`
macro_rules! kff_cast_error {
    ($s:literal) => {
        |_error| {
            error!("{}", $s);
            KFFError::FailedCast
        }
    };
}

// TODO report errors instead of panicking ?
pub fn build_values<FI>(index: &Index<FI>) -> Result<kff::section::Values, KFFError>
where
    FI: FullIndexTrait + Serialize + Sync + Send,
{
    let k = index.get_k();
    let m = index.get_m();
    let ordered = false;
    let max_nb_kmers = (index.get_k() + 1) * 4; // TODO review
    let data_size = std::mem::size_of::<Count>();

    // conversion
    let k = u64::try_from(k).map_err(kff_cast_error!("values of k not fitting into 64 bits"))?;

    let m = u64::try_from(m).map_err(kff_cast_error!("values of m not fitting into 64 bits"))?;
    let ordered = ordered as u64;
    let max_nb_kmers = u64::try_from(max_nb_kmers).map_err(kff_cast_error!(
        "max number of kmers does not fit into 64 bits"
    ))?;
    let data_size =
        u64::try_from(data_size).map_err(kff_cast_error!("usize does not fit into 64 bits"))?;

    let mut values = kff::section::Values::default();
    values.insert("k".to_string(), k);
    values.insert("m".to_string(), m);
    values.insert("ordered".to_string(), ordered);
    values.insert("max".to_string(), max_nb_kmers);
    values.insert("data_size".to_string(), data_size);
    Ok(values)
}

pub fn write_blocks(
    kff_writer: &mut Kff<BufWriter<File>>,
    blocks: &[kff::section::Block],
    m: usize,
    values: &kff::section::Values,
    minimizer: &u64,
) -> Result<(), KFFError> {
    let section = kff::section::Minimizer::new(values).or(Err(KFFError::MissingField))?;
    let bitslice = minimizer.view_bits::<Msb0>();
    let mut bv: BitVec<u8, Msb0> = bitslice.iter().by_vals().collect();

    // TODO please review
    debug_assert!(bv.len() >= 2 * m);
    bv.resize(2 * m, false);

    let bitbox: BitBox<u8, Msb0> = bv.into_boxed_bitslice();

    kff_writer
        .write_minimizer(section, bitbox, blocks)
        .or(Err(KFFError::WriteMinimizer))?;
    Ok(())
}

pub fn create_block(
    context: &str,
    count: &Count,
    minimizer_start_pos: &usize,
    k: usize,
) -> Result<kff::section::Block, KFFError> {
    // build the data vector
    let nb_bases = context.len();
    let nb_kmers = nb_bases - k + 1;
    let data_for_one_kmer = count.to_le_bytes();
    let data_for_all_kmers: Vec<u8> = data_for_one_kmer
        .iter()
        .cycle()
        .take(data_for_one_kmer.len() * nb_kmers)
        .copied()
        .collect();

    let k = u64::try_from(k).map_err(kff_cast_error!("values of k not fitting into 64 bits"))?;
    // build the sequence
    let kmer_seq = kff::Kmer::from_ascii(context.as_bytes(), data_for_all_kmers, ENCODING);
    Ok(kff::section::Block::new(
        k,
        std::mem::size_of::<Count>(),
        kmer_seq,
        *minimizer_start_pos,
    ))
}

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn dump<P: AsRef<Path>, FI: FullIndexTrait + Serialize + Sync + Send>(
    index: &Index<FI>,
    output_file: P,
) -> Result<(), KFFError> {
    index.par_write_kff(output_file)
}
