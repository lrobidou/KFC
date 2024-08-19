use crate::index::FullIndexTrait;
use crate::{Count, Index};
use bitvec::view::BitView;
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use kff::Kff;
use log::error;
use serde::Serialize;
use std::{fs::File, io::BufWriter, path::Path};
use thiserror::Error;

const ENCODING: u8 = 0b00011011;

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
fn build_values<FI>(index: &Index<FI>) -> Result<kff::section::Values, KFFError>
where
    FI: FullIndexTrait + Serialize,
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

fn write_blocks(
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

fn create_block(
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

fn write_all_context<FI>(
    kff_writer: &mut Kff<BufWriter<File>>,
    index: &Index<FI>,
    values: &kff::section::Values,
) -> Result<(), KFFError>
where
    FI: FullIndexTrait + Serialize,
{
    let mut previous_minimizer = None;
    let mut blocks = vec![];
    let k = index.get_k();
    let m = index.get_m();

    for (context, minimizer, minimizer_start_pos, count) in index.context_iterator() {
        match previous_minimizer {
            None => {
                previous_minimizer = Some(minimizer);
            }
            Some(previous_minimizer_val) if previous_minimizer_val == minimizer => {}
            Some(previous_minimizer_val) => {
                // new minimizer, so we write the blocks
                write_blocks(kff_writer, &blocks, m, values, &previous_minimizer_val)?;
                blocks.clear();
                previous_minimizer = Some(minimizer);
            }
        }
        let block = create_block(&context, &count, &minimizer_start_pos, k)?;
        blocks.push(block);
    }

    if !blocks.is_empty() {
        write_blocks(kff_writer, &blocks, m, values, &previous_minimizer.unwrap())?;
    };
    Ok(())
}

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn dump<P: AsRef<Path>, FI: FullIndexTrait + Serialize>(
    index: &Index<FI>,
    output_file: P,
) -> Result<(), KFFError> {
    // TODO discuss: our kmers are non unique and canonical, right ?
    let header = kff::section::Header::new(1, 0, ENCODING, false, true, b"".to_vec())
        .or(Err(KFFError::InvalidHeader))?;
    let mut kff_writer = Kff::create(output_file, header).or(Err(KFFError::CreationFailure))?;

    let values = build_values(index)?;
    kff_writer
        .write_values(values.clone())
        .or(Err(KFFError::WriteValues))?;

    write_all_context(&mut kff_writer, index, &values)?;
    kff_writer.finalize().or(Err(KFFError::FinalizeFailure))?;

    Ok(())
}