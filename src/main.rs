#![cfg_attr(feature = "nightly", feature(core_intrinsics))]
use clap::{Args, Parser, Subcommand};
use serde::bin;
// TODO use better library
use ::serde::Serialize;
use fastxgz::fasta_reads;
use index::{CompleteIndex, FullIndexTrait, Index, StrippedIndex};
use itertools::Itertools;
use mashmap::MashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use superkmers_computation::is_canonical;

type Minimizer = u64;
type HashSuperKmer = u64;
type Count = u16;

mod compute_left_and_right;
mod index;
mod minimizer_iter;
mod serde;
mod superkmer;
mod superkmers_computation;
mod two_bits;

// use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
use superkmer::{reverse_complement, Superkmer};

/// Macro to print the name of a variable and its value
mod macros {
    macro_rules! debug_print {
        ($var:expr) => {{
            let value = $var;
            println!("{} = {:?}", stringify!($var), &value);
            value
        }};
    }

    pub(crate) use debug_print;
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index counting the k-mers of a FASTA/Q file
    Build(BuildArgs),
    Dump(DumpArgs),
    KFFDump(KFFDumpArgs),
}

#[derive(Args, Debug)]
struct BuildArgs {
    /// K-mer size
    #[arg(short)]
    k: usize,
    /// Minimizer size
    #[arg(short, default_value_t = 20)]
    m: usize,
    /// Solidity threshold
    #[arg(short, long, default_value_t = 1)]
    threshold: Count,
    /// Input file (FASTA/Q, possibly gzipped)
    #[arg(short, long)]
    input: String,
    /// Output file (no dump by default)
    #[arg(short, long)]
    output: Option<String>,
    /// Check against the results of KMC
    #[arg(long)]
    check_kmc: Option<String>,
    /// Dump full index, with superkmer informations
    #[arg(short, long)]
    fulldump: bool,
}

#[derive(Args, Debug)]
struct DumpArgs {
    /// Input index
    #[arg(short, long)]
    input_index: String,
    /// Output txt file
    #[arg(long)]
    output_text: Option<String>,
    #[arg(long)]
    output_kff: Option<String>,
}

#[derive(Args, Debug)]
struct KFFDumpArgs {
    /// Input index
    #[arg(short, long)]
    input_kff: String,
    /// Output txt file
    #[arg(short, long)]
    output_text: String,
}

/// Prints some statistics about the index
fn stats<FI: FullIndexTrait + Serialize>(index: &Index<FI>, k: usize) {
    for i in 0..index.get_hyperkmers().get_nb_inserted() {
        let slice = index.get_hyperkmers().get_hyperkmer_from_id(i);
        assert_eq!(slice.len(), k - 1);
    }

    let nb_base_in_large_hyperkmers: usize = index
        .get_large_hyperkmers()
        .iter()
        .map(|large_hk| large_hk.0)
        .sum();
    let number_of_hyperkmers = index.get_hyperkmers().get_nb_inserted();
    let number_of_large_hyperkmers = index.get_large_hyperkmers().len();
    use macros::debug_print as p;
    println!("===== stats =====");
    p!(number_of_hyperkmers);
    p!(number_of_large_hyperkmers);
    println!("nb bases in hyperkmers: {}", number_of_hyperkmers * (k - 1));
    println!(
        "nb base in large hyperkmers: {}",
        nb_base_in_large_hyperkmers
    );
}

/// Query the index with the output of KMC
fn compare_to_kmc<P: AsRef<Path>, FI: FullIndexTrait + Serialize>(
    index: &Index<FI>,
    kmc_file: P,
    k: usize,
    m: usize,
    threshold: Count,
) {
    println!("mass comparison using kmc");
    let file = File::open(kmc_file).unwrap();
    let reader = BufReader::new(file);

    let (mut ok, mut ko) = (0, 0);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split = line.split('\t');
        let kmer = &split.next().unwrap();
        let count = &split.next().unwrap();
        let ground_truth_count: usize = count.parse().unwrap();

        if ground_truth_count >= threshold as usize {
            let kfc_res = index.search_kmer(kmer.as_bytes(), k, m);
            if kfc_res as usize == ground_truth_count {
                ok += 1;
            } else {
                println!("error on {kmer}, {kfc_res} != {ground_truth_count}");
                ko += 1;
            }
        }
    }
    println!("ok = {ok}\nko={ko}");
}

/// Checks it two `MashMap` are equal
fn check_equal_mashmap<K: Ord, V: Ord>(map0: &MashMap<K, V>, map1: &MashMap<K, V>) -> bool {
    // OPTIMIZE this is a naive implementation
    let mut v0 = map0.iter().collect_vec();
    let mut v1 = map1.iter().collect_vec();
    v0.sort();
    v1.sort();
    v0 == v1
}

fn main() {
    let args = Cli::parse();
    simple_logger::SimpleLogger::new().env().init().unwrap();
    match args.command {
        Command::Build(args) => {
            let sequences: Vec<String> = fasta_reads(args.input)
                .unwrap()
                .map(|rcstring| rcstring.to_string())
                .collect();
            let sequences: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();

            let k = args.k;
            let m = args.m;
            let threshold = args.threshold;

            // build the index
            let index = Index::<CompleteIndex>::index(k, m, threshold, &sequences);
            // use KMC's output as a query against the index
            if let Some(kmc_file) = args.check_kmc {
                compare_to_kmc(&index, kmc_file, k, m, threshold);
            }

            stats(&index, k);

            // write the index to disk
            if let Some(output_file) = args.output {
                println!("fulldump : {}", args.fulldump);
                if args.fulldump {
                    serde::bin::dump(&index, &output_file).expect("impossible to dump");

                    #[cfg(debug_assertions)]
                    {
                        let index2 = serde::bin::load(&output_file).expect("impossible to load");
                        debug_assert!(index == index2);
                    }
                } else {
                    // not a full dump: let's remove infos from the index
                    let index = index.remove_superkmer_infos();
                    serde::bin::dump(&index, &output_file).expect("impossible to dump");

                    #[cfg(debug_assertions)]
                    {
                        let index2 = serde::bin::load(&output_file).expect("impossible to load");
                        debug_assert!(index == index2);
                    }
                }
            }
        }
        Command::Dump(args) => {
            let input = args.input_index;
            let index: Index<StrippedIndex> = bin::load(&input).expect("unable to read the index");
            if let Some(kff_path) = args.output_kff {
                serde::kff::dump(&index, kff_path.clone()).unwrap();

                // // TODO remove
                // // read kff and dump to txt
                // let file = std::fs::File::open(kff_path.clone()).expect("unable to read kff file");
                // let buffer = BufReader::new(file);
                // let reader = kff::Kff::read(std::io::BufReader::new(buffer))
                //     .expect("cannot open the kff file");
                // // assert!(reader.check().unwrap());  // surprising behavior, "breaks" the reader
                // println!("{:?}", reader.header());
                // println!("{:?}", reader.values());
                // let kmers = reader.kmers().map(|x| {
                //     println!("{:?}", x);
                //     let kmer = x.unwrap();
                //     let seq = kmer.seq(0b00011110);
                //     let mut data: [u8; 2] = [0, 0];
                //     debug_assert_eq!(kmer.data().len(), 2);
                //     data.copy_from_slice(kmer.data());

                //     (
                //         String::from_utf8(seq).expect("cannot parse utf-8"),
                //         Count::from_le_bytes(data),
                //     )
                // });
                // let file = std::fs::File::create("test_kff.txt").expect("unable to read kff file");
                // let mut buffer = BufWriter::new(file);
                // for (kmer, count) in kmers {
                //     writeln!(buffer, "{}\t{}", kmer, count).unwrap();
                // }

                if let Some(plain_text_path) = args.output_text {
                    serde::plain_text::plain_text(&index, plain_text_path);
                }
            }
        }
        Command::KFFDump(args) => {
            let mut file = kff::Kff::<std::io::BufReader<std::fs::File>>::open(args.input_kff)
                .expect("could not open kff file");
            let output_file = File::create(args.output_text).expect("unable to open output file");
            let mut buffer = BufWriter::with_capacity(9000000, output_file);
            let encoding = *(file.header().encoding());

            while let Some(kmer_section) = file.next_kmer_section() {
                let kmer_section = kmer_section.expect("could not read the kmer section");
                let mut kmers_counts: HashMap<String, u16> = HashMap::new();
                for kmer in kmer_section {
                    // extract count
                    let mut slice = [0; 2];
                    slice.copy_from_slice(kmer.data());
                    let count = Count::from_le_bytes(slice);
                    // compute kmer
                    let kmer = {
                        let seq = kmer.seq(encoding);
                        if is_canonical(&seq) {
                            String::from_utf8(seq).unwrap()
                        } else {
                            reverse_complement(&seq)
                        }
                    };
                    // store kmer
                    kmers_counts
                        .entry(kmer)
                        .and_modify(|x| *x = x.saturating_add(count))
                        .or_insert(count);
                }
                // no more kmer in this section => dump the kmers
                for (kmer, count) in kmers_counts {
                    writeln!(buffer, "{}\t{}", kmer, count).unwrap();
                }
            }
        }
    }
}

// TODO tests
// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_suffix_preffix() {
//         let x = "aazerty______";
//         let y = "aazerty-___---------";
//         assert_eq!(common_prefix_length(x, y), 7);
//         assert_eq!(common_prefix_length(y, x), 7);

//         assert_eq!(common_prefix_length("x", "y"), 0);
//         assert_eq!(common_prefix_length("", ""), 0);
//         assert_eq!(common_prefix_length("xyyyyyyyyyy", "xyy"), 3);

//         // suffix
//         let x = "aazerty_____-erc---------";
//         let y = "aazerty-___erc---------";
//         assert_eq!(common_suffix_length(x, y), 12);
//         assert_eq!(common_suffix_length(y, x), 12);

//         assert_eq!(common_suffix_length("x", "y"), 0);
//         assert_eq!(common_suffix_length("", ""), 0);
//         assert_eq!(common_suffix_length("xyyyyyyyyyy", "xyy"), 2);
//     }

//     #[test]
//     fn test_get_left_and_rigth_extended_hk() {
//         // in the sequence CAGATGGTTCAACCCTTAAGTTAGCGCTTATGGGATCAC
//         //     previous sk :              CCTTAAGTTAGCGCTTATGGGATCAC (GTTAGCGCTT)
//         //     current sk :            ACCCTTAAGTTAGCGCTTATG (CCCTTAAGTT)
//         //        next sk : CAGATGGTTCAACCCTTAAGTTAGCGCTTA (AACCCTTAAG)

//         //    SuperKmerInfos { superkmer: "GTGATCCCATAAGCGCTAACTTAAGG", minimizer: "AAGCGCTAAC", was_read_canonical: false, start_of_minimizer_as_read: 16715, start_of_superkmer_as_read: 16709 }
//         //    SuperKmerInfos { superkmer: "CATAAGCGCTAACTTAAGGGT", minimizer: "AACTTAAGGG", was_read_canonical: false, start_of_minimizer_as_read: 16708, start_of_superkmer_as_read: 16707 }
//         //    SuperKmerInfos { superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA", minimizer: "AACCCTTAAG", was_read_canonical: true, start_of_minimizer_as_read: 16706, start_of_superkmer_as_read: 16696 }
//         //    ("ACTTAAGGGT", "TAAGCGCTAACTTAAGGGT")
//         //    left extended hk: ("AGCGCTAACTTAAGG", 15)
//         //    right extended hk: ("TAAGCGCTAACTTAAGGGT", 10)
//         //    TAAGCGCTAACTTAAGGGTTGAACCATCTG
//         let previous_sk = SuperKmerInfos {
//             superkmer: "GTGATCCCATAAGCGCTAACTTAAGG".into(),
//             minimizer: "AAGCGCTAAC".into(),
//             was_read_canonical: false,
//             start_of_minimizer_as_read: 16715,
//             start_of_superkmer_as_read: 16709,
//         };
//         let current_sk = SuperKmerInfos {
//             superkmer: "CATAAGCGCTAACTTAAGGGT".into(),
//             minimizer: "AACTTAAGGG".into(),
//             was_read_canonical: false,
//             start_of_minimizer_as_read: 16708,
//             start_of_superkmer_as_read: 16707,
//         };
//         let next_sk = SuperKmerInfos {
//             superkmer: "CAGATGGTTCAACCCTTAAGTTAGCGCTTA".into(),
//             minimizer: "AACCCTTAAG".into(),
//             was_read_canonical: true,
//             start_of_minimizer_as_read: 16706,
//             start_of_superkmer_as_read: 16696,
//         };
//         let (l, r) = get_left_and_rigth_extended_hk(&previous_sk, &current_sk, &next_sk);
//         //   left: "TAAGCGCTAACTTAAGGGT"
//         //  right: "ACCCTTAAGTTAGCGCTTA"
//         assert_eq!(r.0, "TAAGCGCTAACTTAAGGGT");
//         assert_eq!(l.0, "CATAAGCGCTAACTTAAGG");
//         //    sk right: ACCCTTAAGTTAGCGCTTA / TAAGCGCTAACTTAAGGGT
//         //    sk left: CCTTAAGTTAGCGCTTATG / CATAAGCGCTAACTTAAGG
//     }

//     #[test]
//     fn test_get_hyperkmer_left_and_rigth_id() {
//         use rand::{distributions::Alphanumeric, Rng}; // 0.8
//         let mut hk_count: HKCount = HKCount::new();
//         let mut hyperkmers: Vec<String> = Vec::new();

//         let nb_insertions = 10000;

//         // random insertions in hyperkmers
//         for _ in 0..nb_insertions {
//             let s: String = rand::thread_rng()
//                 .sample_iter(&Alphanumeric)
//                 .take(100)
//                 .map(char::from)
//                 .collect();
//             hyperkmers.push(get_canonical_kmer(&s).1);
//         }

//         // random minimizers
//         let mut minimizers = Vec::new();
//         for _ in 0..(hyperkmers.len() - 1) {
//             let minimizer: String = rand::thread_rng()
//                 .sample_iter(&Alphanumeric)
//                 .take(20)
//                 .map(char::from)
//                 .collect();
//             minimizers.push(get_canonical_kmer(&minimizer).1);
//         }

//         // link minimizer and hyperkmers
//         // start by inserting random values
//         let mut rng = rand::thread_rng();
//         for _ in 0..3 {
//             for minimizer in &minimizers {
//                 let random_left = rng.gen_range(0..hyperkmers.len());
//                 let random_left_start = rng.gen_range(0..hyperkmers[random_left].len() / 2);
//                 let random_left_end =
//                     rng.gen_range(hyperkmers[random_left].len() / 2..hyperkmers[random_left].len());

//                 let random_right = rng.gen_range(0..hyperkmers.len());
//                 let random_rigth_start = rng.gen_range(0..hyperkmers[random_right].len() / 2);
//                 let random_rigth_end = rng
//                     .gen_range(hyperkmers[random_right].len() / 2..hyperkmers[random_right].len());

//                 let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

//                 hk_count.insert_new_entry_in_hyperkmer_count(
//                     minimizer,
//                     (random_left, random_left_start, random_left_end, true),
//                     (random_right, random_rigth_start, random_rigth_end, true),
//                     random_count,
//                 );
//             }
//         }

//         // then link minimizer and hyperkmers
//         for (i, minimizer) in minimizers.iter().enumerate() {
//             let random_count: u16 = rng.gen_range(0..hyperkmers.len()) as u16;
//             let overlap = minimizer.len();
//             hk_count.insert_new_entry_in_hyperkmer_count(
//                 minimizer,
//                 (i, 0, overlap, true),
//                 (i + 1, 0, overlap, true),
//                 random_count,
//             );
//         }

//         // add another random values
//         for _ in 0..10 {
//             for minimizer in &minimizers {
//                 let random_left = rng.gen_range(0..hyperkmers.len());
//                 let random_left_start = rng.gen_range(0..hyperkmers[random_left].len() / 2);
//                 let random_left_end =
//                     rng.gen_range(hyperkmers[random_left].len() / 2..hyperkmers[random_left].len());

//                 let random_right = rng.gen_range(0..hyperkmers.len());
//                 let random_rigth_start = rng.gen_range(0..hyperkmers[random_right].len() / 2);
//                 let random_rigth_end = rng
//                     .gen_range(hyperkmers[random_right].len() / 2..hyperkmers[random_right].len());

//                 let random_count = rng.gen_range(0..hyperkmers.len()) as u16;

//                 hk_count.insert_new_entry_in_hyperkmer_count(
//                     minimizer,
//                     (random_left, random_left_start, random_left_end, true),
//                     (random_right, random_rigth_start, random_rigth_end, true),
//                     random_count,
//                 );
//             }
//         }

//         // among all the random values, we are still able to get back our real data
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizers[5], &hyperkmers[5]);
//         assert_eq!(hk_id, Some(5));
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizers[5], &hyperkmers[6]);
//         assert_eq!(hk_id, Some(6));

//         // query a non existant minimiser
//         let minimizer: String = rand::thread_rng()
//             .sample_iter(&Alphanumeric)
//             .take(20)
//             .map(char::from)
//             .collect();
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmers[5]);
//         assert_eq!(hk_id, None);
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmers[6]);
//         assert_eq!(hk_id, None);

//         // query an existing minimizer, but with a random hyperkmer
//         let hyperkmer: String = rand::thread_rng()
//             .sample_iter(&Alphanumeric)
//             .take(100)
//             .map(char::from)
//             .collect();
//         let hk_id = hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmer);
//         assert_eq!(hk_id, None);
//         let hk_id = hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmer);
//         assert_eq!(hk_id, None);

//         // query an existing minimizer, but with a random existing hyperkmer (TODO probability of failure is small but not null)
//         let hk_id =
//             hk_count.get_extended_hyperkmer_left_id(&hyperkmers, &minimizer, &hyperkmers[3]);
//         assert_eq!(hk_id, None);
//         let hk_id =
//             hk_count.get_extended_hyperkmer_right_id(&hyperkmers, &minimizer, &hyperkmers[4]);
//         assert_eq!(hk_id, None);
//     }

//     // #[test]
//     // fn test_suf_of_x_is_pref_of_y() {
//     //     assert_eq!(suf_of_x_is_pref_of_y("abcdefghij", "fghijyhy"), 5);
//     // }
// }
