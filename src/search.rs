use crate::{
    compute_left_and_right::get_left_and_rigth_of_sk, extended_hyperkmers::ExtendedHyperkmers,
    superkmers_computation::compute_superkmers_linear_streaming, Count,
};

use super::HKCount;

pub fn search_kmer(
    hk_count: &HKCount,
    hyperkmers: &ExtendedHyperkmers,
    kmer: &[u8],
    k: usize,
    m: usize,
) -> Count {
    let mut sks = compute_superkmers_linear_streaming(kmer, k, m).unwrap();
    let superkmer = &sks.next().unwrap(); // there can be only one
    debug_assert_eq!(sks.next(), None);
    let (left_sk, right_sk) = get_left_and_rigth_of_sk(superkmer);

    hk_count.count_occurence_kmer(
        hyperkmers,
        &superkmer.get_minimizer(),
        &left_sk,
        &right_sk,
        k,
        m,
    )
}

#[cfg(test)]
mod tests {
    use crate::{
        hyperkmers_counts::HKMetadata, superkmer::SubsequenceMetadata, two_bits::encode_minimizer,
    };

    use super::*;

    #[test]
    fn test_search() {
        let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        let mut hk: HKCount = HKCount::new();
        let minimizer = "AAAAAAAAAA";
        let minimizer = encode_minimizer(minimizer.bytes());
        let mut hyperkmers = ExtendedHyperkmers::new(kmer.len(), 7);
        let count = 34;

        let search_result = search_kmer(&hk, &hyperkmers, kmer.as_bytes(), kmer.len(), 10);
        assert!(search_result == 0);
        let s = "TGATGAGTACGTAGCGAAAAAAAAA";
        hyperkmers.add_new_ext_hyperkmer(&SubsequenceMetadata::whole_string(s.as_bytes()));

        let s = "AAAAAAAAAGGGTACGTGCATGCAGTGACGG";
        hyperkmers.add_new_ext_hyperkmer(&SubsequenceMetadata::whole_string(s.as_bytes()));

        hk.insert_new_entry_in_hyperkmer_count(
            &minimizer,
            &HKMetadata {
                index: 0,
                start: 0,
                end: 25,
                change_orientation: false,
            },
            &HKMetadata {
                index: 1,
                start: 0,
                end: 31,
                change_orientation: false,
            },
            count,
        );
        let search_result = search_kmer(&hk, &hyperkmers, kmer.as_bytes(), kmer.len(), 10);
        assert!(search_result == count);
    }
}
