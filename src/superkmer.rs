use crate::two_bits::decode_2bits;
use crate::Minimizer;

use super::superkmers_computation::is_canonical;
use super::two_bits;
use itertools::Itertools;
use std::iter::Map;
use std::iter::Rev;
use xxhash_rust::const_xxh3::xxh3_64;

const REVCOMP_TAB_CHAR: [char; 255] = {
    let mut tab = ['A'; 255];
    tab[b'A' as usize] = 'T';
    tab[b'T' as usize] = 'A';
    tab[b'C' as usize] = 'G';
    tab[b'G' as usize] = 'C';
    tab
};

pub fn reverse_complement(seq: &[u8]) -> String {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB_CHAR.get_unchecked(*base as usize) })
        .collect()
}

// TODO duplication: move to a module
const REVCOMP_TAB: [u8; 255] = {
    let mut tab = [0; 255];
    tab[b'A' as usize] = b'T';
    tab[b'T' as usize] = b'A';
    tab[b'C' as usize] = b'G';
    tab[b'G' as usize] = b'C';
    tab
};

// TODO "discuss" is there no copy here ? What is the cost of moving references ?
pub fn reverse_complement_no_copy(
    seq: impl DoubleEndedIterator<Item = u8>,
) -> Map<Rev<impl DoubleEndedIterator<Item = u8>>, fn(u8) -> u8> {
    seq.rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(base as usize) })
}

// states of SubsequenceMetadata
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BitPacked {
    total_base_in_sequence: usize,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoBitPacked;

// OPTIMIZE duplication of `read` when used in Superkmer
/// Represents a subsequence, possibly in reverse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubsequenceMetadata<'a, Packing> {
    read: &'a [u8],
    start: usize,
    end: usize,
    same_orientation: bool,
    packing: Packing,
}

impl<'a> SubsequenceMetadata<'a, NoBitPacked> {
    pub fn new(read: &'a [u8], start: usize, end: usize, same_orientation: bool) -> Self {
        debug_assert!(start <= read.len());
        debug_assert!(end <= read.len());
        Self {
            read,
            start, // in base
            end,   // in base
            same_orientation,
            packing: NoBitPacked {},
        }
    }

    // pub fn whole_string(read: &'a [u8]) -> Self {
    //     Self::new(read, 0, read.len(), true)
    // }

    pub fn is_canonical(&self) -> bool {
        let subsequence = &self.read[self.start..self.end];
        let is_original_subsequence_canonical = is_canonical(subsequence);
        self.same_orientation == is_original_subsequence_canonical
    }

    pub fn to_canonical(self) -> Self {
        let subsequence = &self.read[self.start..self.end];

        if is_canonical(subsequence) == self.same_orientation {
            self
        } else {
            Self {
                read: self.read,
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        }
    }

    pub fn to_canonical_string(&self) -> String {
        let subsequence = &self.read[self.start..self.end];
        if is_canonical(subsequence) {
            String::from_utf8(Vec::from(subsequence)).unwrap()
        } else {
            reverse_complement(subsequence)
        }
    }

    // pub fn equal_str(&self, other: &str) -> bool {
    //     self.to_string() == other
    // }

    // pub fn equal(&self, other: &SubsequenceMetadata<NoBitPacked>) -> bool {
    //     if self.len() != other.len() {
    //         false
    //     } else {
    //         self.common_prefix_length(other) == other.len()
    //     }
    // }

    // pub fn ends_with(&self, other: &SubsequenceMetadata<NoBitPacked>) -> bool {
    //     if self.len() < other.len() {
    //         false
    //     } else {
    //         self.common_suffix_length(other) == other.len()
    //     }
    // }

    // pub fn starts_with(&self, other: &SubsequenceMetadata<NoBitPacked>) -> bool {
    //     if self.len() < other.len() {
    //         false
    //     } else {
    //         self.common_prefix_length(other) == other.len()
    //     }
    // }

    // pub fn encode_as_minimizer(&self) -> Vec<u64> {
    //     let subsequence = &self.read[self.start..self.end];
    //     if self.same_orientation {
    //         two_bits::encode_minimizer(subsequence.iter().copied())
    //     } else {
    //         two_bits::encode_2bits_u64(
    //             reverse_complement_no_copy(subsequence.iter().copied()),
    //             self.len(),
    //         )
    //     }
    // }

    pub fn dump_as_2bits(&self, slice: &mut [u8]) {
        let subsequence = &self.read[self.start..self.end];
        if self.same_orientation {
            let iter = two_bits::encode_2bits(subsequence.iter().copied(), self.len());

            for (ret, src) in slice.iter_mut().zip(iter) {
                *ret = src;
            }
        } else {
            let iter = two_bits::encode_2bits(
                reverse_complement_no_copy(subsequence.iter().copied()),
                self.len(),
            );
            for (ret, src) in slice.iter_mut().zip(iter) {
                *ret = src;
            }
        }
    }

    // pub fn common_prefix_length(&self, other: &SubsequenceMetadata<NoBitPacked>) -> usize {
    //     if self.same_orientation && other.same_orientation {
    //         let mut x = self.read[self.start..self.end].iter().copied();
    //         let mut y = other.read[other.start..other.end].iter().copied();
    //         iter_prefix_len(x, y)
    //     } else if self.same_orientation && !other.same_orientation {
    //         let mut x = self.read[self.start..self.end].iter().copied();
    //         let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
    //         iter_prefix_len(&mut x, &mut y)
    //     } else if !self.same_orientation && other.same_orientation {
    //         let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
    //         let mut y = other.read[other.start..other.end].iter().copied();
    //         iter_prefix_len(&mut x, &mut y)
    //     } else {
    //         // !self.same_orientation && !other.same_orientation
    //         let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
    //         let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
    //         iter_prefix_len(&mut x, &mut y)
    //     }
    // }

    // pub fn common_suffix_length(&self, other: &SubsequenceMetadata<NoBitPacked>) -> usize {
    //     if self.same_orientation && other.same_orientation {
    //         let mut x = self.read[self.start..self.end].iter().copied();
    //         let mut y = other.read[other.start..other.end].iter().copied();
    //         iter_suffix_len(&mut x, &mut y)
    //     } else if self.same_orientation && !other.same_orientation {
    //         let mut x = self.read[self.start..self.end].iter().copied();
    //         let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
    //         iter_suffix_len(&mut x, &mut y)
    //     } else if !self.same_orientation && other.same_orientation {
    //         let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
    //         let mut y = other.read[other.start..other.end].iter().copied();
    //         iter_suffix_len(&mut x, &mut y)
    //     } else {
    //         //f !self.same_orientation && !other.same_orientation {
    //         let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
    //         let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
    //         iter_suffix_len(&mut x, &mut y)
    //     }
    // }

    pub fn common_prefix_length_with_bitpacked(
        &self,
        other: &SubsequenceMetadata<BitPacked>,
    ) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.read[self.start..self.end].iter().copied();
            let mut y = other.decode_2bits();
            iter_prefix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.read[self.start..self.end].iter().copied();
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_prefix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let self_iterator = self.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let mut y = other.decode_2bits();
            iter_prefix_len(&mut x, &mut y)
        } else {
            // !self.same_orientation && !other.same_orientation
            let self_iterator = self.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_prefix_len(&mut x, &mut y)
        }
    }

    pub fn common_suffix_length_with_bitpacked(
        &self,
        other: &SubsequenceMetadata<BitPacked>,
    ) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.read[self.start..self.end].iter().copied();
            let mut y = other.decode_2bits();
            iter_suffix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.read[self.start..self.end].iter().copied();
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_suffix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let self_iterator = self.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let mut y = other.decode_2bits();
            iter_suffix_len(&mut x, &mut y)
        } else {
            // !self.same_orientation && !other.same_orientation
            let self_iterator = self.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_suffix_len(&mut x, &mut y)
        }
    }

    pub fn equal_bitpacked(&self, other: &SubsequenceMetadata<BitPacked>) -> bool {
        if self.len() != other.len() {
            false
        } else {
            self.common_prefix_length_with_bitpacked(other) == other.len()
        }
    }

    pub fn hash(&self) -> u64 {
        let subsequence = &self.read[self.start..self.end];
        if self.same_orientation {
            xxh3_64(subsequence)
        } else {
            let revcomp = reverse_complement_no_copy(subsequence.iter().copied());
            // TODO copy into a vec
            xxh3_64(&revcomp.collect_vec())
        }
    }
}

impl<'a> std::fmt::Display for SubsequenceMetadata<'a, NoBitPacked> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let subsequence = &self.read[self.start..self.end];
        let string = if self.same_orientation {
            String::from_utf8(Vec::from(subsequence)).unwrap()
        } else {
            reverse_complement(subsequence)
        };
        write!(f, "{}", string)
    }
}

impl<'a> SubsequenceMetadata<'a, BitPacked> {
    pub fn whole_bitpacked(bytes: &'a [u8], nb_bases: usize) -> Self {
        debug_assert!((nb_bases / 4) + ((nb_bases % 4 != 0) as usize) == bytes.len());
        Self {
            read: bytes,
            start: 0,
            end: nb_bases,
            same_orientation: true,
            packing: BitPacked {
                total_base_in_sequence: nb_bases,
            },
        }
    }

    pub fn ends_with_nobitpacked(&self, other: &SubsequenceMetadata<NoBitPacked>) -> bool {
        if self.len() < other.len() {
            false
        } else {
            other.common_suffix_length_with_bitpacked(self) == other.len()
        }
    }

    pub fn starts_with_nobitpacked(&self, other: &SubsequenceMetadata<NoBitPacked>) -> bool {
        if self.len() < other.len() {
            false
        } else {
            other.common_prefix_length_with_bitpacked(self) == other.len()
        }
    }

    pub fn decode_2bits(&self) -> impl DoubleEndedIterator<Item = u8> + 'a {
        debug_assert!(self.end <= self.packing.total_base_in_sequence);
        debug_assert!(self.read.len() * 4 >= self.packing.total_base_in_sequence);
        #[cfg(debug_assertions)]
        {
            // test to check that the reverse is working
            // TODO do more check and tests for the reverse decoding
            let truth = decode_2bits(
                self.read.iter().cloned(),
                self.start,
                self.end,
                self.packing.total_base_in_sequence,
            )
            .collect_vec();
            let what_i_made = decode_2bits(
                self.read.iter().cloned(),
                self.start,
                self.end,
                self.packing.total_base_in_sequence,
            )
            .rev()
            .collect_vec()
            .iter()
            .rev()
            .copied()
            .collect_vec();

            debug_assert_eq!(truth.len(), self.end - self.start);

            debug_assert_eq!(truth, what_i_made);
        }
        decode_2bits(
            self.read.iter().cloned(),
            self.start,
            self.end,
            self.packing.total_base_in_sequence,
        )
    }
}

impl<'a> std::fmt::Display for SubsequenceMetadata<'a, BitPacked> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let iter = decode_2bits(
            self.read.iter().copied(),
            self.start,
            self.end,
            self.packing.total_base_in_sequence,
        );
        let s = String::from_utf8(iter.collect_vec()).expect("Invalid string");
        write!(f, "{}", s)
    }
}

impl<'a, Packing> SubsequenceMetadata<'a, Packing>
where
    Packing: Copy,
{
    pub fn change_orientation(&self) -> Self {
        Self {
            read: self.read,
            start: self.start,
            end: self.end,
            same_orientation: !self.same_orientation,
            packing: self.packing,
        }
    }

    pub fn change_orientation_if(&self, cond: bool) -> Self {
        if cond {
            Self {
                read: self.read,
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        } else {
            *self
        }
    }

    pub fn start(&self) -> usize {
        self.start
    }
    pub fn end(&self) -> usize {
        self.end
    }
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Extract a subsequence
    /// When not bit packed, equivalent to:
    /// ```
    /// Self::whole_string(self.to_string()[start..end])
    /// ```
    pub fn subsequence(&self, start: usize, end: usize) -> SubsequenceMetadata<Packing> {
        if self.same_orientation {
            Self {
                read: self.read,
                start: self.start + start,
                end: self.start + end,
                same_orientation: self.same_orientation,
                packing: self.packing,
            }
        } else {
            Self {
                read: self.read,
                start: self.end - end,
                end: self.end - start,
                same_orientation: self.same_orientation,
                packing: self.packing,
            }
        }
    }
}

fn iter_prefix_len(mut x: impl Iterator<Item = u8>, mut y: impl Iterator<Item = u8>) -> usize {
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}

fn iter_suffix_len(
    x: &mut impl DoubleEndedIterator<Item = u8>,
    y: &mut impl DoubleEndedIterator<Item = u8>,
) -> usize {
    let mut x = x.rev();
    let mut y = y.rev();
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}
#[derive(Debug, Clone, Copy, PartialEq)]

pub struct Superkmer<'a> {
    pub read: &'a [u8],
    minimizer: u64,
    start_mini: usize,
    end_mini: usize,
    pub superkmer: SubsequenceMetadata<'a, NoBitPacked>,
}

impl<'a> Superkmer<'a> {
    pub fn new(
        read: &'a [u8],
        start_mini: usize,
        end_mini: usize,
        start_sk: usize,
        end_sk: usize,
        same_orientation: bool,
    ) -> Self {
        let minizer_subsequence = &read[start_mini..end_mini];
        let minimizer = if same_orientation {
            two_bits::encode_minimizer(minizer_subsequence.iter().copied())
        } else {
            two_bits::encode_minimizer(reverse_complement_no_copy(
                minizer_subsequence.iter().copied(),
            ))
        };

        Self {
            read,
            minimizer,
            start_mini,
            end_mini,
            superkmer: SubsequenceMetadata::<NoBitPacked>::new(
                read,
                start_sk,
                end_sk,
                same_orientation,
            ),
        }
    }

    // pub fn new_with_hash(
    //     read: &'a [u8],
    //     start_mini: usize,
    //     end_mini: usize,
    //     start_sk: usize,
    //     end_sk: usize,
    //     hash_minimizer: u64,
    //     same_orientation: bool,
    // ) -> Self {
    //     let minizer_subsequence = &read[start_mini..end_mini];
    //     let minimizer = if same_orientation {
    //         two_bits::encode_minimizer(minizer_subsequence.iter().copied())
    //     } else {
    //         two_bits::encode_minimizer(reverse_complement_no_copy(
    //             minizer_subsequence.iter().copied(),
    //         ))
    //     };

    //     Self {
    //         read,
    //         minimizer,
    //         start_mini,
    //         end_mini,
    //         superkmer: SubsequenceMetadata::<NoBitPacked>::new(
    //             read,
    //             start_sk,
    //             end_sk,
    //             same_orientation,
    //         ),
    //     }
    // }

    pub fn hash_superkmer(&self) -> u64 {
        self.superkmer.hash()
    }

    pub fn get_minimizer(&self) -> Minimizer {
        self.minimizer
    }

    pub fn start_of_minimizer(&self) -> usize {
        self.start_mini
    }

    pub fn end_of_minimizer(&self) -> usize {
        self.end_mini
    }

    // pub fn m(&self) -> usize {
    //     self.minimizer.len()
    // }

    pub fn is_canonical_in_the_read(&self) -> bool {
        self.superkmer.same_orientation
    }

    // pub fn print(&self) -> String {
    //     format!(
    //         "{} ({} {})",
    //         String::from_utf8(Vec::from(
    //             &self.superkmer.read[self.superkmer.start..self.superkmer.end]
    //         ))
    //         .unwrap(),
    //         self.minimizer.to_string(),
    //         self.minimizer.same_orientation
    //     )
    // }
}

mod tests {
    use two_bits::encode_2bits;

    use super::*;

    #[test]
    pub fn test_dump() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let hk = SubsequenceMetadata::new(bytes, 0, bytes.len(), true);

        let dest = &mut [0, 0, 0];
        hk.to_canonical().dump_as_2bits(dest);

        let expected: Vec<u8> = encode_2bits(bytes.iter().copied(), read.len()).collect();
        assert_eq!(dest, expected.as_slice())
    }
}
