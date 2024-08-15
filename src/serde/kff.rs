use crate::{Count, Index};
use bitvec::view::BitView;
use bitvec::{boxed::BitBox, order::Msb0, vec::BitVec};
use itertools::Itertools;
use kff::kmer::seq2bits;
use kff::Kff;
use log::error;
use thiserror::Error;

use std::{fs::File, io::BufWriter, path::Path};

type KffWriter = kff::Kff<std::io::BufWriter<std::fs::File>>;

const ENCODING: u8 = 0b00011011;
/// Creates a kff from a header
fn header<P: AsRef<Path>>(output_file: P) -> KffWriter {
    // TODO discuss: our kmers are non unique and canonical, right ?
    let header = kff::section::Header::new(1, 0, ENCODING, false, true, b"".to_vec())
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
fn build_values(index: &Index) -> Result<kff::section::Values, KFFError> {
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
    // let mut bv: BitVec<u8, Msb0> = bv.into_vec().iter().take(bytes2store(m)).collect();
    debug_assert!(bv.len() >= 2 * m);
    bv.resize(2 * m, false);

    let bitbox: BitBox<u8, Msb0> = bv.into_boxed_bitslice();
    // println!(
    //     "writing minimizer {:?}",
    //     String::from_utf8(kff::kmer::bits2seq(&bitbox, ENCODING)).unwrap()
    // );

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

fn write_all_context(
    kff_writer: &mut Kff<BufWriter<File>>,
    index: &Index,
    values: &kff::section::Values,
) -> Result<(), KFFError> {
    let mut previous_minimizer = None;
    let mut blocks = vec![];
    let k = index.get_k();
    let m = index.get_m();

    for (context, minimizer, minimizer_start_pos, count) in index.context_iterator() {
        // println!("context{}", context);
        // println!(
        //     "minimizer = {}",
        //     &context[minimizer_start_pos..minimizer_start_pos + m]
        // );
        // println!("minimizer_start_pos = {minimizer_start_pos}");
        // println!("count = {count}");

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

fn write_large_context(
    kff_writer: &mut Kff<BufWriter<File>>,
    index: &Index,
    values: &kff::section::Values,
) -> Result<(), KFFError> {
    let mut previous_minimizer = None;
    let mut blocks = vec![];
    let k = index.get_k();
    let m = index.get_m();

    for (context, minimizer, minimizer_start_pos, count) in index.large_context_iterator() {
        match previous_minimizer {
            None => {
                previous_minimizer = Some(minimizer);
            }
            Some(previous_minimizer_val) if previous_minimizer_val == minimizer => {}
            Some(_) => {
                // new minimizer, so we write the blocks
                write_blocks(kff_writer, &blocks, m, values, &minimizer)?;
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
pub fn dump<P: AsRef<Path>>(index: &Index, output_file: P) -> Result<(), KFFError> {
    let mut kff_writer = header(output_file);

    let values_normal_context = build_values(index)?;
    kff_writer
        .write_values(values_normal_context.clone())
        .or(Err(KFFError::WriteValues))?;

    write_all_context(&mut kff_writer, index, &values_normal_context)?;
    kff_writer.finalize().expect("impossible to finalize");

    Ok(())
}

// fn ceil_to_8(n: usize) -> usize {
//     (n + 7) & !(7)
// }

// fn bytes2store(k: usize) -> usize {
//     ceil_to_8(k * 2) / 8
// }
// Dump all kmers and their abundance from the index into `output_file`.
// K-mers are not sorted.
// pub fn buggy_dump<P: AsRef<Path>>(index: &Index, output_file: P) -> Result<(), KFFError> {
//     // encoding: [b'A', b'C', b'T', b'G']
//     let header = kff::section::Header::new(1, 0, ENCODING, false, true, b"".to_vec())
//         .expect("invalid header");
//     let mut kff_writer = Kff::create(output_file, header).expect("unable to initiate kff");
//     let values_normal_context = build_values_normal_context(index)?;
//     kff_writer
//         .write_values(values_normal_context.clone())
//         .or(Err(KFFError::WriteValues))?;

//     let contexts = index.normal_context_iterator().take(1).collect_vec();
//     let (context, minimizer, start_minimizer_pos, count) = contexts[0].clone();

//     let nb_bases = context.len();
//     let nb_kmers = nb_bases - index.get_k() + 1;
//     let data_for_one_kmer = count.to_le_bytes();
//     let mut data_for_all_kmers: Vec<u8> = data_for_one_kmer
//         .iter()
//         .cycle()
//         .take(data_for_one_kmer.len() * nb_kmers)
//         .copied()
//         .collect();
//     let size = data_for_all_kmers.len();
//     data_for_all_kmers[0] = 0;
//     data_for_all_kmers[2] = 2;
//     data_for_all_kmers[5] = 5;
//     data_for_all_kmers[6] = 6;
//     data_for_all_kmers[9] = 9;
//     data_for_all_kmers[size - 1] = 100;

//     // build the sequence
//     let kmer_seq = kff::Kmer::from_ascii(context.as_bytes(), data_for_all_kmers.clone(), ENCODING);
//     let kmer_seq2 = kff::Kmer::from_ascii(context.as_bytes(), data_for_all_kmers, ENCODING);

//     let k = index.get_k() as u64;

//     // push the sequence into the blocks
//     let block = kff::section::Block::new(
//         k,
//         std::mem::size_of::<Count>(),
//         kmer_seq,
//         start_minimizer_pos,
//     );

//     let block2 = kff::section::Block::new(
//         k,
//         std::mem::size_of::<Count>(),
//         kmer_seq2,
//         start_minimizer_pos,
//     );
//     let blocks = vec![block, block2];

//     let section =
//         kff::section::Minimizer::new(&values_normal_context).or(Err(KFFError::MissingField))?;

//     let bitslice = minimizer.view_bits::<Msb0>();
//     let bv: BitVec<u8, Msb0> = bitslice.iter().by_vals().collect();
//     let mut bv: BitVec<u8, Msb0> = bv
//         .into_vec()
//         .iter()
//         .take(bytes2store(index.get_m()))
//         .collect();
//     bv.resize(2 * index.get_m(), false);
//     let bitbox: BitBox<u8, Msb0> = bv.into_boxed_bitslice();
//     kff_writer
//         .write_minimizer(section, bitbox, &blocks)
//         .or(Err(KFFError::WriteMinimizer))?;
//     kff_writer.finalize().expect("impossible to finalize");

//     Ok(())
// }
