#[macro_use]
pub mod components;
mod computation;
mod extraction;

use super::Minimizer;
use crate::serde::kff::{build_values, create_block, write_blocks, KFFError, ENCODING};
use crate::Count;
use crate::{
    compute_left_and_right::get_left_and_rigth_of_sk,
    superkmers_computation::compute_superkmers_linear_streaming,
};
use extraction::{extract_context, extract_kmers_from_contexts_associated_to_a_minimizer};

use components::{HKCount, LargeExtendedHyperkmer, ParallelExtendedHyperkmers, SuperKmerCounts};

use crate::buckets::Buckets;
use computation::{first_stage, second_stage};
use kff::Kff;
use rayon::prelude::*;
use serde::{
    de::{SeqAccess, Visitor},
    ser::SerializeStruct,
    Deserialize, Deserializer, Serialize, Serializer,
};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::RwLock;
use std::sync::{Arc, Mutex};
use std::time::Instant;
// -- type state pattern
pub trait FullIndexTrait: std::marker::Sync {}

#[derive(Serialize, Deserialize, PartialEq)]
pub struct CompleteIndex {
    super_kmer_counts: Buckets<SuperKmerCounts>,
    discarded_minimizers: Buckets<HashMap<Minimizer, u16>>,
}

#[derive(Serialize, Deserialize, PartialEq)]
pub struct StrippedIndex {}

impl FullIndexTrait for CompleteIndex {}
impl FullIndexTrait for StrippedIndex {}

/// Index of KFC
/// The index holds the hyperkmers and their counts.
/// During contruction, the index also holds the counts of superkmers and information about minimizers.
pub struct Index<FI: FullIndexTrait + Sync + Send + Serialize> {
    hk_count: Buckets<HKCount>,
    hyperkmers: Arc<RwLock<ParallelExtendedHyperkmers>>,
    /// vector of larger extended hyperkmers // TODO document
    large_hyperkmers: Arc<RwLock<Vec<LargeExtendedHyperkmer>>>,
    k: usize,
    m: usize,
    superkmers_infos: FI,
}

impl<FI: FullIndexTrait + Sync + Send + Serialize> Serialize for Index<FI> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let hyperkmers = self.hyperkmers.read().unwrap();
        let large_hyperkmers = self.large_hyperkmers.read().unwrap();
        let k = self.k;
        let m = self.m;
        // let superkmers_infos = self.superkmers_infos;

        let mut state = serializer.serialize_struct("Index", 6)?;
        state.serialize_field("hk_count", &self.hk_count)?;
        state.serialize_field("hyperkmers", &*hyperkmers)?;
        state.serialize_field("large_hyperkmers", &*large_hyperkmers)?;
        state.serialize_field("k", &k)?;
        state.serialize_field("m", &m)?;
        state.serialize_field("superkmers_infos", &self.superkmers_infos)?;
        // serialize other fields
        state.end()
    }
}

impl<'a, FI: FullIndexTrait + Sync + Send + Serialize + for<'de> Deserialize<'de>> Deserialize<'a>
    for Index<FI>
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'a>,
    {
        deserializer.deserialize_struct(
            "Index",
            &[
                "hk_count",
                "hyperkmers",
                "large_hyperkmers",
                "k",
                "m",
                "superkmers_infos",
            ],
            IndexVisitor::<FI> {
                _phantom_data: PhantomData::<FI>,
            },
        )
    }
}

struct IndexVisitor<FI> {
    _phantom_data: PhantomData<FI>,
}

impl<'a, FI: FullIndexTrait + Sync + Send + Serialize + for<'de> Deserialize<'de>> Visitor<'a>
    for IndexVisitor<FI>
{
    type Value = Index<FI>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct Index")
    }

    fn visit_seq<V>(self, mut seq: V) -> Result<Self::Value, V::Error>
    where
        V: SeqAccess<'a>,
    {
        // TODO is the len arg correct?
        let hk_count = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
        let hyperkmers = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(1, &self))?;
        let large_hyperkmers = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(2, &self))?;
        let k = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(3, &self))?;
        let m = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(4, &self))?;
        let superkmers_infos = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(5, &self))?;

        // let hk_count = Arc::new(RwLock::new(hk_count));
        let hyperkmers = Arc::new(RwLock::new(hyperkmers));
        let large_hyperkmers = Arc::new(RwLock::new(large_hyperkmers));

        // let superkmers_infos = Arc::new(RwLock::new(superkmers_infos));
        Ok(Index {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            k,
            m,
            superkmers_infos,
        })
    }
}

macro_rules! eq_arc_rwlock {
    ($x:expr, $y:expr) => {{
        *$x.read().expect("could acquire read lock") == *$y.read().expect("could acquire read lock")
    }};
}

impl<FI: FullIndexTrait + PartialEq + Send + Sync + Serialize + Serialize> PartialEq for Index<FI> {
    fn eq(&self, other: &Self) -> bool {
        self.k == other.k
            && self.m == other.m
            && self.hk_count == other.hk_count
            && eq_arc_rwlock!(self.hyperkmers, other.hyperkmers)
            && eq_arc_rwlock!(self.large_hyperkmers, other.large_hyperkmers)
            && self.superkmers_infos == other.superkmers_infos
    }
}

impl Index<CompleteIndex> {
    /// Constructs a new `Index` indexing a set of sequences
    #[allow(clippy::self_named_constructors)] // Self named constructor ? I want it that way 🎵
    pub fn index<P: AsRef<Path>>(k: usize, m: usize, threshold: Count, path: P) -> Self {
        let start_fisrt_stage = Instant::now();
        let (super_kmer_counts, hk_count, hyperkmers, large_hyperkmers) =
            first_stage(&path, k, m, threshold);
        let time_first_stage = start_fisrt_stage.elapsed().as_secs();
        println!(
            "time first stage: {} second{}",
            time_first_stage,
            if time_first_stage > 1 { "s" } else { "" }
        );

        let start_second_stage = Instant::now();
        let discarded_minimizers = second_stage(
            &super_kmer_counts,
            &hk_count,
            hyperkmers.clone(),
            large_hyperkmers.clone(),
            &path,
            k,
            m,
            threshold,
        );
        let time_second_stage = start_second_stage.elapsed().as_secs();
        println!(
            "time second stage: {} second{}",
            time_second_stage,
            if time_second_stage > 1 { "s" } else { "" }
        );

        Self {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            superkmers_infos: CompleteIndex {
                super_kmer_counts,
                discarded_minimizers,
            },
            k,
            m,
        }
    }

    pub fn remove_superkmer_infos(self) -> Index<StrippedIndex> {
        Index {
            hk_count: self.hk_count,
            hyperkmers: self.hyperkmers,
            large_hyperkmers: self.large_hyperkmers,
            k: self.k,
            m: self.m,
            superkmers_infos: StrippedIndex {},
        }
    }
}

impl<FI> Index<FI>
where
    FI: FullIndexTrait + Serialize + Sync + Send + Serialize,
{
    #[cfg(test)]
    pub fn new(
        hk_count: Buckets<HKCount>,
        hyperkmers: ParallelExtendedHyperkmers,
        large_hyperkmers: Vec<LargeExtendedHyperkmer>,
        superkmers_infos: FI,
        k: usize,
        m: usize,
    ) -> Self {
        let hyperkmers = Arc::new(RwLock::new(hyperkmers));
        let large_hyperkmers = Arc::new(RwLock::new(large_hyperkmers));
        Self {
            hk_count,
            hyperkmers,
            large_hyperkmers,
            superkmers_infos,
            k,
            m,
        }
    }

    pub fn get_hyperkmers(&self) -> &Arc<RwLock<ParallelExtendedHyperkmers>> {
        &self.hyperkmers
    }

    pub fn get_large_hyperkmers(&self) -> &Arc<RwLock<Vec<LargeExtendedHyperkmer>>> {
        &self.large_hyperkmers
    }

    pub fn search_kmer(&self, kmer: &[u8], k: usize, m: usize) -> Count {
        // there can be only one superkmer
        let mut sks = compute_superkmers_linear_streaming(kmer, k, m).unwrap();
        let superkmer = &sks.next().unwrap();
        debug_assert_eq!(sks.next(), None);

        let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

        let hk_count_for_minimizer = self.hk_count.get_from_id_u64(superkmer.get_minimizer());

        let hk_count_for_minimizer = hk_count_for_minimizer
            .read()
            .expect("could not acquire read lock");

        let hyperkmers = self.hyperkmers.read().expect("could not acquire read lock");
        let large_hyperkmers = self
            .large_hyperkmers
            .read()
            .expect("could not acquire read lock");

        hk_count_for_minimizer.count_occurence_kmer(
            &hyperkmers,
            &large_hyperkmers,
            &superkmer.get_minimizer(),
            &left_sk,
            &right_sk,
            k,
            m,
        )
    }

    pub fn par_write_kmers<P: AsRef<Path>>(&self, kmer_threshold: Count, path: P) {
        // create file and wrap it in mutex
        let file = File::create(path).unwrap();
        let buffer = BufWriter::with_capacity(100_000, file);
        let buffer = Arc::new(Mutex::new(buffer));

        let hk_count_chunks = self.hk_count.chunks();

        hk_count_chunks.par_iter().for_each(|chunk| {
            let k = self.k;
            let m = self.m;
            let hyperkmers = self.hyperkmers.read().unwrap();
            let large_hyperkmers = self.large_hyperkmers.read().unwrap();
            let hk_count_chunk = chunk.read().unwrap();

            let mut km_counts_grouped_by_key =
                hk_count_chunk.get_data().iter_group_by_key().peekable();

            // OPTIMIZE: possibility to change this allocation
            let mut kmers_and_count = Vec::with_capacity(hk_count_chunk.get_data().len());
            while let Some(kmers) = extract_kmers_from_contexts_associated_to_a_minimizer(
                &mut km_counts_grouped_by_key,
                &hyperkmers,
                &large_hyperkmers,
                &k,
                &m,
            ) {
                for kmer_and_count in kmers {
                    kmers_and_count.push(kmer_and_count);
                }
            }

            let mut buffer = buffer.as_ref().lock().unwrap();
            for (kmer, count) in kmers_and_count {
                if count >= kmer_threshold {
                    buffer
                        .write_all(&kmer)
                        .expect("could not write to the file");
                    writeln!(buffer, "\t{}", count).unwrap();
                }
            }
        });
    }

    pub fn par_write_kff<P: AsRef<Path>>(&self, path: P) -> Result<(), KFFError> {
        // TODO discuss: our kmers are non unique and canonical, right ?
        let header = kff::section::Header::new(1, 0, ENCODING, false, true, b"".to_vec())
            .or(Err(KFFError::InvalidHeader))?;
        let mut kff_writer = Kff::create(path, header).or(Err(KFFError::CreationFailure))?;

        let values = build_values(self)?;
        kff_writer
            .write_values(values.clone())
            .or(Err(KFFError::WriteValues))?;

        let kff_writer = Arc::new(Mutex::new(kff_writer));

        let hk_count_chunks = self.hk_count.chunks();

        hk_count_chunks.par_iter().for_each(|chunk| {
            let k = self.k;
            let m = self.m;
            let hyperkmers = self.hyperkmers.read().unwrap();
            let large_hyperkmers = self.large_hyperkmers.read().unwrap();
            let hk_count_chunk = chunk.read().unwrap();

            let minimizers = hk_count_chunk.minimizer_set();

            for minimizer in minimizers {
                let mut blocks = vec![];
                let hk_entries = hk_count_chunk.get_data().get_iter(&minimizer);
                for entry in hk_entries {
                    let (context, minimizer_start_pos) =
                        extract_context(entry, m, &hyperkmers, &large_hyperkmers);
                    let kff_block =
                        create_block(&context, &entry.2, &minimizer_start_pos, k).unwrap();
                    blocks.push(kff_block)
                }
                // TODO a bit stupid to write now
                let mut kff_writer = kff_writer.lock().unwrap();
                write_blocks(&mut kff_writer, &blocks, m, &values, &minimizer).unwrap();
            }
        });
        // TODO do we have finished
        let mut kff_writer = kff_writer.lock().unwrap();
        kff_writer.finalize().or(Err(KFFError::FinalizeFailure))?;
        Ok(())
    }

    pub fn get_k(&self) -> usize {
        self.k
    }

    pub fn get_m(&self) -> usize {
        self.m
    }
}

#[cfg(test)]
mod tests {

    use crate::serde::bin;

    use crate::{
        compute_left_and_right, subsequence::Subsequence,
        superkmers_computation::compute_superkmers_linear_streaming,
    };

    use super::*;
    use components::HKMetadata;
    use itertools::Itertools;

    #[test]
    fn test_search_empty_index() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 10;

        // nothing inserted => nothing is found
        let empty_index: Index<CompleteIndex> = Index::new(
            Buckets::<HKCount>::new(HKCount::new),
            ParallelExtendedHyperkmers::new(k, 5),
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new),
                discarded_minimizers: Buckets::<HashMap<u64, u16>>::new(HashMap::new),
            },
            k,
            m,
        );
        let search_result = empty_index.search_kmer(kmer.as_bytes(), k, m);
        assert!(search_result == 0);
    }

    #[test]
    fn test_search() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 10;

        let superkmers = compute_superkmers_linear_streaming(kmer.as_bytes(), kmer.len(), m)
            .unwrap()
            .collect_vec();
        assert!(superkmers.len() == 1);
        let superkmer = superkmers[0];

        // let's assume the minimizer is as follow
        assert_eq!(superkmer.start_of_minimizer(), 24);
        assert_eq!(superkmer.end_of_minimizer(), 34);

        let hyperkmers = ParallelExtendedHyperkmers::new(kmer.len(), 7);
        let large_hyperkmers = Vec::new();
        let count = 34;

        // computing the left and right context to insert them in vector of hyperkmer
        let (left, right) = compute_left_and_right::get_left_and_rigth_of_sk(&superkmer);
        let size_left = left.len();
        let size_right = right.len();

        assert_eq!(size_left, superkmer.end_of_minimizer() - 1);
        assert_eq!(
            size_right,
            superkmer.superkmer.len() - (superkmer.start_of_minimizer() + 1)
        );

        // adding 'A' to complete the hyperkmers
        let mut left = left.to_string();
        while left.len() < k - 1 {
            left.push('A');
        }
        let mut right = right.to_string();
        while right.len() < k - 1 {
            right.push('A');
        }

        // computing the left and right extended hyperkmers
        let left = Subsequence::new(left.as_bytes(), 0, left.len(), true);
        let right = Subsequence::new(right.as_bytes(), 0, right.len(), true);

        // inserting the hyperkmers
        let (bucket_left, index_left) = hyperkmers.add_new_ext_hyperkmer(&left);
        let (bucket_right, index_right) = hyperkmers.add_new_ext_hyperkmer(&right);

        let left_hk_metadata = HKMetadata::new(
            bucket_left,
            index_left,
            0,
            size_left,
            false,
            left.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let right_hk_metadata = HKMetadata::new(
            bucket_right,
            index_right,
            0,
            size_right,
            false,
            right.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let hk_count = Buckets::<HKCount>::new(HKCount::new);
        let hk_count_chunk = hk_count.get_from_id_u64(superkmer.get_minimizer());
        let mut hk_count_chunk = hk_count_chunk.write().unwrap();
        hk_count_chunk.insert_new_entry_in_hyperkmer_count(
            &superkmer.get_minimizer(),
            &left_hk_metadata,
            &right_hk_metadata,
            count,
        );
        drop(hk_count_chunk);

        let index: Index<CompleteIndex> = Index::new(
            hk_count,
            hyperkmers,
            large_hyperkmers,
            CompleteIndex {
                super_kmer_counts: Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new),
                discarded_minimizers: Buckets::<HashMap<u64, u16>>::new(HashMap::new),
            },
            k,
            m,
        );
        let search_result = index.search_kmer(kmer.as_bytes(), kmer.len(), 10);
        assert_eq!(search_result, count);
    }

    #[test]
    fn test_empty_serde() {
        let filename = "test_empty_serde";
        let k = 10;
        let m = 3;

        let empty_index: Index<CompleteIndex> = Index::new(
            Buckets::<HKCount>::new(HKCount::new),
            ParallelExtendedHyperkmers::new(k, 5),
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new),
                discarded_minimizers: Buckets::<HashMap<u64, u16>>::new(HashMap::new),
            },
            k,
            m,
        );
        let empty_index = empty_index.remove_superkmer_infos();
        bin::dump(&empty_index, filename).unwrap();
        assert!(empty_index == bin::load(filename).unwrap());
    }

    #[test]
    fn test_serde() {
        let filename = "test_serde";
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let k = kmer.len();
        let m = 3;

        let superkmers = compute_superkmers_linear_streaming(kmer.as_bytes(), kmer.len(), 10)
            .unwrap()
            .collect_vec();
        assert!(superkmers.len() == 1);
        let superkmer = superkmers[0];

        // let's assume the minimizer is as follow
        assert_eq!(superkmer.start_of_minimizer(), 24);
        assert_eq!(superkmer.end_of_minimizer(), 34);

        let hyperkmers = ParallelExtendedHyperkmers::new(kmer.len(), 7);
        let count = 34;

        // computing the left and right context to insert them in vector of hyperkmer
        let (left, right) = compute_left_and_right::get_left_and_rigth_of_sk(&superkmer);
        let size_left = left.len();
        let size_right = right.len();

        assert_eq!(size_left, superkmer.end_of_minimizer() - 1);
        assert_eq!(
            size_right,
            superkmer.superkmer.len() - (superkmer.start_of_minimizer() + 1)
        );

        // adding 'A' to complete the hyperkmers
        let mut left = left.to_string();
        while left.len() < k - 1 {
            left.push('A');
        }
        let mut right = right.to_string();
        while right.len() < k - 1 {
            right.push('A');
        }

        // computing the left and right extended hyperkmers
        let left = Subsequence::new(left.as_bytes(), 0, left.len(), true);
        let right = Subsequence::new(right.as_bytes(), 0, right.len(), true);

        // inserting the hyperkmers
        let (bucket_left, index_left) = hyperkmers.add_new_ext_hyperkmer(&left);
        let (bucket_right, index_right) = hyperkmers.add_new_ext_hyperkmer(&right);

        let left_hk_metadata = HKMetadata::new(
            bucket_left,
            index_left,
            0,
            size_left,
            false,
            left.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let right_hk_metadata = HKMetadata::new(
            bucket_right,
            index_right,
            0,
            size_right,
            false,
            right.is_canonical() != superkmer.is_canonical_in_the_read(),
        );

        let hk_count = Buckets::<HKCount>::new(HKCount::new);
        let hk_count_chunk = hk_count.get_from_id_u64(superkmer.get_minimizer());
        let mut hk_count_chunk = hk_count_chunk.write().unwrap();
        hk_count_chunk.insert_new_entry_in_hyperkmer_count(
            &superkmer.get_minimizer(),
            &left_hk_metadata,
            &right_hk_metadata,
            count,
        );
        drop(hk_count_chunk);

        let index: Index<CompleteIndex> = Index::new(
            hk_count,
            hyperkmers,
            Vec::new(),
            CompleteIndex {
                super_kmer_counts: Buckets::<SuperKmerCounts>::new(SuperKmerCounts::new),
                discarded_minimizers: Buckets::<HashMap<u64, u16>>::new(HashMap::new),
            },
            k,
            m,
        );

        let index = index.remove_superkmer_infos();
        bin::dump(&index, filename).unwrap();
        assert!(index == bin::load(filename).unwrap());
    }

    // #[test]
    // fn test_iter_kmers_empty() {
    //     let index: Index<CompleteIndex> = Index::<CompleteIndex>::index(21, 11, 1, &vec![]);
    //     let kmers_iter = index.iter_kmers().flatten();
    //     let kmers = kmers_iter.collect_vec();
    //     assert_eq!(kmers, vec![]);
    // }

    // #[test]
    // fn test_iter_kmers_one_element() {
    //     let v: HashSet<(String, Count)> = HashSet::from_iter(
    //         vec![
    //             ("TAGCTTCTCGCTATTAGCTTC", 1),
    //             ("AGCTTCTCGCTATTAGCTTCA", 1),
    //             ("GCTTCTCGCTATTAGCTTCAA", 1),
    //             ("CTTCTCGCTATTAGCTTCAAT", 1),
    //             ("TTCTCGCTATTAGCTTCAATG", 1),
    //             ("TCTCGCTATTAGCTTCAATGA", 1),
    //             ("CTCGCTATTAGCTTCAATGAT", 1),
    //             ("TCGCTATTAGCTTCAATGATA", 1),
    //             ("CGCTATTAGCTTCAATGATAC", 1),
    //             ("GCTATTAGCTTCAATGATACG", 1),
    //         ]
    //         .into_iter()
    //         .map(|x| (x.0.into(), x.1)),
    //     );

    //     let index = Index::index(21, 11, 1, &vec!["TAGCTTCTCGCTATTAGCTTCAATGATACG"]);
    //     let kmers_results: HashSet<(String, Count)> = HashSet::from_iter(index.iter_kmers());
    //     assert_eq!(kmers_results, v);
    // }
}
