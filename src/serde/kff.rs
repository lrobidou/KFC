use crate::{Count, Index};
use bitvec::view::BitView;
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use itertools::Itertools;
use kff::{Kff, Kmer};
use log::error;
use thiserror::Error;

use std::{error, fs::File, io::BufWriter, path::Path};

type KffWriter = kff::Kff<std::io::BufWriter<std::fs::File>>;

/// Creates a kff from a header
fn header<P: AsRef<Path>>(output_file: P) -> KffWriter {
    // TODO discuss: our kmers are non unique and canonical, right ?
    // let header = section::Header::new(1, 0, 0b00011011, true, true, )?;

    let header = kff::section::Header::new(1, 0, 0b00011110, false, true, b"".to_vec())
        .expect("invalid header");
    Kff::create(output_file, header).expect("unable to initiate kff")
}

macro_rules! kff_cast_error {
    ($s:literal) => {
        |_error| {
            error!("{}", $s);
            KFFError::FailedCast
        }
    };
}

// TODO report errors instead of panicking ?
// we set max to k-1, so that we can store normal contexts here.
fn build_values_normal_context(index: &Index) -> Result<kff::section::Values, KFFError> {
    let k = index.get_k();
    let m = index.get_m();
    let ordered = false;
    let max_nb_kmers = index.get_k() - index.get_m() + 1;
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

#[derive(Error, Debug)]
pub enum KFFError {
    #[error("Failed to cast values to a KFF-compatible format")]
    FailedCast,
    #[error("Failed to write values on a KFF file")]
    WriteValues, // #[error("the data for key `{0}` is not available")]
    // Redaction(String),
    // #[error("invalid header (expected {expected:?}, found {found:?})")]
    // InvalidHeader { expected: String, found: String },
    // #[error("unknown data store error")]
    // Unknown,
    #[error("Missing field")]
    MissingField,
    #[error("Could not write the minimizer section")]
    WriteMinimizer,
}

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn dump<P: AsRef<Path>>(index: &Index, output_file: P) -> Result<(), KFFError> {
    // encoding: [b'A', b'C', b'T', b'G']
    let mut kff_writer = header(output_file);

    let values_normal_context = build_values_normal_context(index)?;
    kff_writer
        .write_values(values_normal_context.clone())
        .or(Err(KFFError::WriteValues))?;

    // kff_writer.write_minimizer(section, minimizer, blocks)

    // let minimizer_iter = index.iter_minimizers();
    // let section =
    //     kff::section::Minimizer::new(&values_normal_context).or(Err(KFFError::MissingField))?;

    let mut previous_minimizer = None;
    let mut blocks = vec![];
    for (context, minimizer, mini_pos, count) in index.normal_context_iterator() {
        match previous_minimizer {
            None => previous_minimizer = Some(minimizer),
            Some(previous_minimizer) => {
                if previous_minimizer == minimizer {
                    // build the data vector
                    let nb_bases = context.len() + index.get_m();
                    let nb_kmers = nb_bases - index.get_k() + 1;
                    let data_for_one_kmer = count.to_le_bytes();
                    let data_for_all_kmers = data_for_one_kmer
                        .iter()
                        .cycle()
                        .take(data_for_one_kmer.len() * nb_kmers)
                        .copied()
                        .collect();

                    // build the sequence
                    let kmer_seq =
                        kff::Kmer::from_ascii(context.as_bytes(), data_for_all_kmers, 0b00011110);

                    let k = u64::try_from(index.get_k())
                        .map_err(kff_cast_error!("values of k not fitting into 64 bits"))?;

                    // push the sequence into the blocks
                    blocks.push(kff::section::Block::new(
                        k,
                        std::mem::size_of::<Count>(),
                        kmer_seq,
                        mini_pos,
                    ));
                } else {
                    println!("dumping");
                    let section = kff::section::Minimizer::new(&values_normal_context)
                        .or(Err(KFFError::MissingField))?;
                    let bitslice = previous_minimizer.view_bits::<Msb0>();
                    let bv: BitVec<u8, Msb0> = bitslice.iter().by_vals().collect();
                    let bitbox: BitBox<u8, Msb0> = bv.into_boxed_bitslice();
                    kff_writer
                        .write_minimizer(section, bitbox, &blocks)
                        .or(Err(KFFError::WriteMinimizer))?;
                    kff_writer.finalize().expect("impossible to finalize");
                    // return Ok(());
                    blocks.clear();
                }
            }
        }
    }

    Ok(())
}

/// Dump all kmers and their abundance from the index into `output_file`.
/// K-mers are not sorted.
pub fn buggy_dump<P: AsRef<Path>>(index: &Index, output_file: P) -> Result<(), KFFError> {
    // encoding: [b'A', b'C', b'T', b'G']
    let header = kff::section::Header::new(1, 0, 0b00011110, false, true, b"".to_vec())
        .expect("invalid header");
    let mut kff_writer = Kff::create(output_file, header).expect("unable to initiate kff");
    let values_normal_context = build_values_normal_context(index)?;
    kff_writer
        .write_values(values_normal_context.clone())
        .or(Err(KFFError::WriteValues))?;

    let contexts = index.normal_context_iterator().take(1).collect_vec();
    let (context, minimizer, start_minimizer_pos, count) = contexts[0].clone();

    let nb_bases = context.len() + index.get_m();
    let nb_kmers = nb_bases - index.get_k() + 1;
    let data_for_one_kmer = count.to_le_bytes();
    let data_for_all_kmers = data_for_one_kmer
        .iter()
        .cycle()
        .take(data_for_one_kmer.len() * nb_kmers)
        .copied()
        .collect();

    // build the sequence
    let kmer_seq = kff::Kmer::from_ascii(context.as_bytes(), data_for_all_kmers, 0b00011110);

    println!("context {:?}", context);

    let k = index.get_k() as u64;

    // push the sequence into the blocks
    let blocks = vec![kff::section::Block::new(
        k,
        std::mem::size_of::<Count>(),
        kmer_seq,
        start_minimizer_pos,
    )];

    let section =
        kff::section::Minimizer::new(&values_normal_context).or(Err(KFFError::MissingField))?;
    let bitslice = minimizer.view_bits::<Msb0>();
    let bv: BitVec<u8, Msb0> = bitslice.iter().by_vals().collect();
    let bitbox: BitBox<u8, Msb0> = bv.into_boxed_bitslice();
    kff_writer
        .write_minimizer(section, bitbox, &blocks)
        .or(Err(KFFError::WriteMinimizer))?;
    kff_writer.finalize().expect("impossible to finalize");

    Ok(())
}
