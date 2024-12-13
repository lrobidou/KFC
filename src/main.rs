#![cfg_attr(feature = "nightly", feature(core_intrinsics))]

use ::serde::Serialize;
use clap::{Args, Parser, Subcommand};
use index::{CompleteIndex, FullIndexTrait, Index, StrippedIndex};
use itertools::Itertools;
use macros::p;
use mashmap::MashMap;
use serde::bin;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use superkmers_computation::is_canonical;

type Minimizer = u64;
type HashSuperKmer = u64;
type Count = u16;

mod buckets;
mod codec;
mod complexity;
mod compute_left_and_right;
mod index;
mod minimizer_iter;
mod ringbuf;
mod serde;
mod simd;
mod subsequence;
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

    pub(crate) use debug_print as p;
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
    /// Dump an index into a text or KFF file
    Dump(DumpArgs),
    /// Dump a KFF file into a text file (warning: KFC only handles KFF files produced by KFC)
    KFFDump(KFFDumpArgs),
}

#[derive(Args, Debug)]
struct BuildArgs {
    /// K-mer size
    #[arg(short)]
    k: usize,
    /// Minimizer size
    #[arg(short, default_value_t = 21)]
    m: usize,
    /// Expert parameter: solidity threshold for superkmers
    #[arg(short, long, default_value_t = 2)]
    threshold_superkmer: Count,
    /// Input file (FASTA/Q, possibly gzipped)
    #[arg(short, long)]
    input: String,
    /// Output file ({input}.kff by default)
    #[arg(short, long)]
    output: Option<String>,
    /// Number of threads (use all core by default)
    #[arg(short = 'T', long)]
    threads: Option<usize>,
    /// Check against the results of KMC (no check by default)
    #[arg(long)]
    check_kmc: Option<String>,
    /// Expert parameter: dump additional information with the index. KFC cannot read these indexes from the CLI. (partial dump by default)
    #[arg(short, long)]
    fulldump: bool,
    /// Print some statistics about the index (no print by default)
    #[arg(short, long)]
    print_stats: bool,
}

#[derive(Args, Debug)]
struct DumpArgs {
    /// Input index
    #[arg(short, long)]
    input_index: String,
    /// Output txt file (no txt output by default)
    #[arg(long)]
    output_text: Option<String>,
    /// Output kff file (no kff output by default)
    #[arg(long)]
    output_kff: Option<String>,
    /// Minimum abundance of the k-mers to output (if text output)
    #[arg(short = 't', long)]
    kmer_threshold: Count,
    /// Number of threads (use all core by default)
    #[arg(short = 'T', long)]
    threads: Option<usize>,
}

#[derive(Args, Debug)]
struct KFFDumpArgs {
    /// Input index
    #[arg(short, long)]
    input_kff: String,
    /// Output txt file
    #[arg(short, long)]
    output_text: String,
    /// Minimum abundance of the k-mers to output
    #[arg(short = 't', long)]
    kmer_threshold: Count,
    /// Number of threads (use all core by default)
    #[arg(short = 'T', long)]
    threads: Option<usize>,
}

/// Formats a u64 in a String
fn format_number_u64(number: u64) -> String {
    number
        .to_string()
        .as_bytes()
        .rchunks(3)
        .rev()
        .map(std::str::from_utf8)
        .collect::<Result<Vec<&str>, _>>()
        .unwrap()
        .join(",")
}

/// Formats a usize in a String
fn format_number_usize(number: usize) -> String {
    number
        .to_string()
        .as_bytes()
        .rchunks(3)
        .rev()
        .map(std::str::from_utf8)
        .collect::<Result<Vec<&str>, _>>()
        .unwrap()
        .join(",")
}

/// Checks a file exists. Exits the program if the path is not a file.
fn check_file_exists<P>(filepath: P)
where
    P: AsRef<Path>,
{
    if filepath.as_ref().exists() {
        if !filepath.as_ref().is_file() {
            let filepath = filepath.as_ref().to_str().unwrap();
            eprintln!("Error: '{filepath}' is not a file.");
            std::process::exit(1);
        }
    } else {
        let filepath = filepath.as_ref().to_str().unwrap();
        eprintln!("Error: file '{filepath}' does not exist.");
        std::process::exit(1);
    }
}

/// Prints some statistics about the index
fn print_stats<FI: FullIndexTrait + Serialize + Sync + Send + Serialize>(
    index: &Index<FI>,
    k: usize,
) {
    let hyperkmers = index.get_hyperkmers();
    let hyperkmers = hyperkmers.read().expect("could not acquire read lock");

    let large_hyperkmers = index.get_large_hyperkmers();
    let large_hyperkmers = large_hyperkmers
        .read()
        .expect("could not acquire read lock");

    let number_of_hyperkmers = hyperkmers.get_nb_inserted();
    let number_of_large_hyperkmers = large_hyperkmers.len();
    let nb_base_in_hyperkmers: usize = number_of_hyperkmers * (k - 1);
    let nb_base_in_large_hyperkmers: usize =
        large_hyperkmers.iter().map(|large_hk| large_hk.0).sum();
    let (nb_minimizers, nb_kmers) = index.count_minimizers_and_kmers();

    println!("===== stats =====");
    println!(
        "indexed {} kmers with {} minimizers",
        format_number_usize(nb_minimizers),
        format_number_usize(nb_kmers)
    );
    println!(
        "number of hyperkmers: {}",
        format_number_usize(number_of_hyperkmers + number_of_large_hyperkmers)
    );

    #[cfg(debug_assertions)]
    println!(
        "    [debug] including {} large hyperkmers",
        format_number_usize(number_of_large_hyperkmers)
    );

    println!(
        "number of bases in hyperkmers: {}",
        format_number_usize(nb_base_in_hyperkmers + nb_base_in_large_hyperkmers)
    );
    #[cfg(debug_assertions)]
    println!(
        "    [debug] including {} bases in large hyperkmers",
        format_number_usize(nb_base_in_large_hyperkmers)
    );

    #[cfg(debug_assertions)]
    {
        println!(
            "    [debug] minimizers =  {:?}",
            index.get_nb_minimizers_in_each_bucket()
        );
        println!(
            "    [debug] entries =  {:?}",
            index.get_nb_entries_in_each_bucket()
        );
    }
}

/// Query the index with the output of KMC
fn compare_to_kmc<P: AsRef<Path>, FI: FullIndexTrait + Serialize + Sync + Send + Serialize>(
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

// TODO implement this in mashmap ?
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
    // simple_logger::SimpleLogger::new()
    //     .with_level(log::LevelFilter::Warn)
    //     .env()
    //     .init()
    //     .unwrap();
    match args.command {
        Command::Build(args) => {
            let k = args.k;
            let m = args.m;
            let threshold = args.threshold_superkmer;

            // check that k and m are odd
            if k % 2 != 1 {
                eprintln!("Error: k must be odd.");
                std::process::exit(1);
            }
            if m % 2 != 1 {
                eprintln!("Error: m must be odd.");
                std::process::exit(1);
            }

            check_file_exists(&args.input);
            if let Some(ref kmc_file) = args.check_kmc {
                check_file_exists(kmc_file);
            }

            // set the number of threads
            if let Some(nb_threads) = args.threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(nb_threads)
                    .build_global()
                    .expect("unable to configure the threads");
            }

            // get output filename
            let output_file = match args.output {
                Some(filename) => filename,
                None => {
                    let mut path = PathBuf::from(&args.input);
                    // Set the new extension or add one if there's none
                    path.set_extension("kfc");
                    path.to_string_lossy().to_string()
                }
            };

            // build the index
            let index = Index::<CompleteIndex>::index(k, m, threshold, args.input);
            // use KMC's output as a query against the index
            if let Some(kmc_file) = args.check_kmc {
                compare_to_kmc(&index, kmc_file, k, m, threshold);
            }

            if args.print_stats {
                print_stats(&index, k);
            }

            // write the index to disk
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
        Command::Dump(args) => {
            let input = args.input_index;
            let kmer_threshold = args.kmer_threshold;

            check_file_exists(&input);

            // set the number of threads
            if let Some(nb_threads) = args.threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(nb_threads)
                    .build_global()
                    .expect("unable to configure the threads");
            }

            let index: Index<StrippedIndex> = bin::load(&input).expect("unable to read the index");
            if let Some(kff_path) = args.output_kff {
                serde::kff::dump(&index, kff_path.clone()).unwrap();
            }
            if let Some(plain_text_path) = args.output_text {
                serde::plain_text::plain_text(&index, kmer_threshold, plain_text_path);
            }
        }
        Command::KFFDump(args) => {
            let kmer_threshold = args.kmer_threshold;

            check_file_exists(&args.input_kff);

            // set the number of threads
            if let Some(nb_threads) = args.threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(nb_threads)
                    .build_global()
                    .expect("unable to configure the threads");
            }

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
                    if count >= kmer_threshold {
                        writeln!(buffer, "{}\t{}", kmer, count).unwrap();
                    }
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
