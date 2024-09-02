use itertools::Itertools;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::{
    superkmer::{
        reverse_complement, reverse_complement_ascii_to_ascii, reverse_complement_no_copy,
    },
    two_bits::{self, decode_2bits},
};

use super::superkmers_computation::is_canonical;

// states of Subsequence
#[derive(Debug, Clone, PartialEq)]
pub struct BitPacked<'a> {
    total_base_in_sequence: usize,
    read: &'a [u8],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoBitPacked<'a> {
    read: &'a [u8],
}

// OPTIMIZE duplication of `read` when used in Superkmer
/// Represents a subsequence, possibly in reverse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subsequence<Packing> {
    start: usize,
    end: usize,
    same_orientation: bool,
    packing: Packing,
}

impl<'a> Subsequence<NoBitPacked<'a>> {
    pub fn new(read: &'a [u8], start: usize, end: usize, same_orientation: bool) -> Self {
        debug_assert!(start <= read.len());
        debug_assert!(end <= read.len());
        debug_assert!(start <= end);
        Self {
            start, // in base
            end,   // in base
            same_orientation,
            packing: NoBitPacked { read },
        }
    }

    // pub fn whole_string(read: &'a [u8]) -> Self {
    //     Self::new(read, 0, read.len(), true)
    // }

    pub fn is_canonical(&self) -> bool {
        let subsequence = &self.packing.read[self.start..self.end];
        let is_original_subsequence_canonical = is_canonical(subsequence);
        self.same_orientation == is_original_subsequence_canonical
    }

    pub fn to_canonical(self) -> Self {
        let subsequence = &self.packing.read[self.start..self.end];

        if is_canonical(subsequence) == self.same_orientation {
            self
        } else {
            Self {
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        }
    }

    pub fn dump_as_2bits(&self, slice: &mut [u8]) {
        let subsequence = &self.packing.read[self.start..self.end];
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

    pub fn common_prefix_length_with_bitpacked(&self, other: &Subsequence<BitPacked>) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.packing.read[self.start..self.end].iter().copied();
            let mut y = other.decode_2bits();
            iter_prefix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.packing.read[self.start..self.end].iter().copied();
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_prefix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let self_iterator = self.packing.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let mut y = other.decode_2bits();
            iter_prefix_len(&mut x, &mut y)
        } else {
            // !self.same_orientation && !other.same_orientation
            let self_iterator = self.packing.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_prefix_len(&mut x, &mut y)
        }
    }

    pub fn common_suffix_length_with_bitpacked(&self, other: &Subsequence<BitPacked>) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.packing.read[self.start..self.end].iter().copied();
            let mut y = other.decode_2bits();
            iter_suffix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.packing.read[self.start..self.end].iter().copied();
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_suffix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let self_iterator = self.packing.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let mut y = other.decode_2bits();
            iter_suffix_len(&mut x, &mut y)
        } else {
            // !self.same_orientation && !other.same_orientation
            let self_iterator = self.packing.read[self.start..self.end].iter().copied();
            let mut x = reverse_complement_no_copy(self_iterator);
            let other_iterator = other.decode_2bits();
            let mut y = reverse_complement_no_copy(other_iterator);
            iter_suffix_len(&mut x, &mut y)
        }
    }

    pub fn equal_bitpacked(&self, other: &Subsequence<BitPacked>) -> bool {
        if self.len() != other.len() {
            false
        } else {
            self.common_prefix_length_with_bitpacked(other) == other.len()
        }
    }

    pub fn hash(&self) -> u64 {
        let subsequence = &self.packing.read[self.start..self.end];
        if self.same_orientation {
            xxh3_64(subsequence)
        } else {
            let revcomp = reverse_complement_no_copy(subsequence.iter().copied());
            // TODO copy into a vec
            xxh3_64(&revcomp.collect_vec())
        }
    }
}

impl<'a> std::fmt::Display for Subsequence<NoBitPacked<'a>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let subsequence = &self.packing.read[self.start..self.end];
        let string = if self.same_orientation {
            String::from_utf8(Vec::from(subsequence)).unwrap()
        } else {
            reverse_complement(subsequence)
        };
        write!(f, "{}", string)
    }
}

impl<'a> Subsequence<NoBitPacked<'a>> {
    pub fn as_vec(&self) -> Vec<u8> {
        let subsequence = &self.packing.read[self.start..self.end];
        if self.same_orientation {
            Vec::from(subsequence)
        } else {
            reverse_complement_ascii_to_ascii(subsequence)
        }
    }
}

impl<'a> Subsequence<BitPacked<'a>> {
    pub fn whole_bitpacked(read: &'a [u8], nb_bases: usize) -> Self {
        #[cfg(debug_assertions)]
        {
            // let read = lock.get_slice_from_internal_id(id);
            debug_assert!((nb_bases / 4) + ((nb_bases % 4 != 0) as usize) == read.len());
        }
        Self {
            start: 0,
            end: nb_bases,
            same_orientation: true,
            packing: BitPacked {
                total_base_in_sequence: nb_bases,
                read,
            },
        }
    }

    pub fn ends_with_nobitpacked(&self, other: &Subsequence<NoBitPacked>) -> bool {
        if self.len() < other.len() {
            false
        } else {
            other.common_suffix_length_with_bitpacked(self) == other.len()
        }
    }

    pub fn starts_with_nobitpacked(&self, other: &Subsequence<NoBitPacked>) -> bool {
        if self.len() < other.len() {
            false
        } else {
            other.common_prefix_length_with_bitpacked(self) == other.len()
        }
    }

    pub fn decode_2bits(&self) -> impl DoubleEndedIterator<Item = u8> + '_ {
        debug_assert!(self.end <= self.packing.total_base_in_sequence);
        debug_assert!(self.packing.read.len() * 4 >= self.packing.total_base_in_sequence);
        #[cfg(debug_assertions)]
        {
            // test to check that the reverse is working
            // TODO do more check and tests for the reverse decoding
            let truth = decode_2bits(
                self.packing.read.iter().cloned(),
                self.start,
                self.end,
                self.packing.total_base_in_sequence,
            )
            .collect_vec();
            let what_i_made = decode_2bits(
                self.packing.read.iter().cloned(),
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
            self.packing.read.iter().cloned(),
            self.start,
            self.end,
            self.packing.total_base_in_sequence,
        )
    }
}

impl<'a> std::fmt::Display for Subsequence<BitPacked<'a>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let iter = decode_2bits(
            self.packing.read.iter().copied(),
            self.start,
            self.end,
            self.packing.total_base_in_sequence,
        )
        .collect_vec();
        let string = if self.same_orientation {
            String::from_utf8(iter).unwrap()
        } else {
            reverse_complement(&iter)
        };
        write!(f, "{}", string)
    }
}

// TODO no copy (?)
impl<'a> Subsequence<BitPacked<'a>> {
    pub fn as_vec(&self) -> Vec<u8> {
        let iter = decode_2bits(
            self.packing.read.iter().copied(),
            self.start,
            self.end,
            self.packing.total_base_in_sequence,
        )
        .collect_vec();
        if self.same_orientation {
            iter
        } else {
            reverse_complement_ascii_to_ascii(&iter)
        }
    }
}

impl<Packing> Subsequence<Packing>
// where
// Packing: Copy,
{
    pub fn change_orientation(self) -> Self {
        Self {
            start: self.start,
            end: self.end,
            same_orientation: !self.same_orientation,
            packing: self.packing,
        }
    }

    pub fn change_orientation_if(self, cond: bool) -> Self {
        if cond {
            Self {
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        } else {
            self
        }
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn same_orientation(&self) -> bool {
        self.same_orientation
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Extract a subsequence
    /// When not bit packed, equivalent to:
    /// ```
    /// Self::whole_string(self.to_string()[start..end])
    /// ```
    pub fn subsequence(self, start: usize, end: usize) -> Subsequence<Packing> {
        if self.same_orientation {
            Self {
                start: self.start + start,
                end: self.start + end,
                same_orientation: self.same_orientation,
                packing: self.packing,
            }
        } else {
            Self {
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
