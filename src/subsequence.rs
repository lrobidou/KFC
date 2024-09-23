use itertools::Itertools;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::{
    codec::{
        align_left_iterator::AlignLeftIterator,
        prefix::{
            common_prefix_length_ff, common_prefix_length_fr, common_prefix_length_rf,
            common_prefix_length_rr, RevCompEncoder,
        },
        suffix::{
            common_suffix_length_ff, common_suffix_length_fr, common_suffix_length_rf,
            common_suffix_length_rr,
        },
        Decoder, Encoder,
    },
    superkmer::{
        reverse_complement, reverse_complement_ascii_to_ascii, reverse_complement_no_copy,
    },
};

use super::superkmers_computation::is_canonical;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

// states of Subsequence
#[derive(Debug, Clone, PartialEq)]
pub struct BitPacked<'a> {
    total_base_in_sequence: usize,
    read: &'a [u64],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoBitPacked<'a> {
    read: &'a [u8],
}

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

    pub fn get_read(&self) -> &[u8] {
        self.packing.read
    }

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

    pub fn dump_as_2bits(&self, slice: &mut [u64]) {
        let subsequence = &self.packing.read[self.start..self.end];
        if self.same_orientation {
            let encoder = Encoder::new(subsequence);
            for (ret, src) in slice.iter_mut().zip(encoder) {
                *ret = src;
            }
        } else {
            let encoder = RevCompEncoder::new(subsequence);
            for (ret, src) in slice.iter_mut().zip(encoder) {
                *ret = src;
            }
        }
    }

    pub fn common_prefix_length_with_bitpacked(&self, other: &Subsequence<BitPacked>) -> usize {
        let self_read = &self.packing.read[self.start..self.end];
        let other_read = other.packing.read;
        let other_read = AlignLeftIterator::new(other_read, other.start, other.end).collect_vec();

        let common_prefix_len = if self.same_orientation && other.same_orientation {
            common_prefix_length_ff(self_read, &other_read, other.len())
        } else if self.same_orientation && !other.same_orientation {
            common_prefix_length_fr(self_read, &other_read, other.len())
        } else if !self.same_orientation && other.same_orientation {
            common_prefix_length_rf(self_read, &other_read, other.len())
        } else {
            common_prefix_length_rr(self_read, &other_read, other.len())
        };

        // check that I made no programming mistake
        #[cfg(debug_assertions)]
        {
            let prefix = prefix_str(&self.to_string(), &other.to_string());
            debug_assert_eq!(prefix, common_prefix_len);
        }
        common_prefix_len
    }

    pub fn common_suffix_length_with_bitpacked(&self, other: &Subsequence<BitPacked>) -> usize {
        let self_read = &self.packing.read[self.start..self.end];
        let other_read = other.packing.read;
        let other_read = AlignLeftIterator::new(other_read, other.start, other.end).collect_vec();

        let common_suffix_len = if self.same_orientation && other.same_orientation {
            common_suffix_length_ff(self_read, &other_read, other.len())
        } else if self.same_orientation && !other.same_orientation {
            common_suffix_length_fr(self_read, &other_read, other.len())
        } else if !self.same_orientation && other.same_orientation {
            common_suffix_length_rf(self_read, &other_read, other.len())
        } else {
            common_suffix_length_rr(self_read, &other_read, other.len())
        };

        #[cfg(debug_assertions)]
        {
            let suffix = suffix_str(&self.to_string(), &other.to_string());
            debug_assert_eq!(suffix, common_suffix_len);
        }
        common_suffix_len
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
    pub fn whole_bitpacked(read: &'a [u64], nb_bases: usize) -> Self {
        #[cfg(debug_assertions)]
        {
            // let read = lock.get_slice_from_internal_id(id);
            debug_assert!((nb_bases / 32) + ((nb_bases % 32 != 0) as usize) == read.len());
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

    // pub fn decode_2bits(&self) -> Vec<u8> {
    //     let decoder = Decoder::new(self.packing.read, self.packing.total_base_in_sequence);
    //     let ascii_whole = decoder.collect_vec();
    //     let ascii = ascii_whole[self.start..self.end]
    //         .iter()
    //         .copied()
    //         .collect_vec();
    //     ascii
    // }
}

impl<'a> std::fmt::Display for Subsequence<BitPacked<'a>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let ascii = self.as_vec();
        let string = String::from_utf8(ascii).unwrap();
        write!(f, "{}", string)
    }
}

impl<'a> Subsequence<BitPacked<'a>> {
    pub fn as_vec(&self) -> Vec<u8> {
        // TODO remove when the Decoder can handle different starts
        // TODO we decode more than necessary
        let bytes = Decoder::new(self.packing.read, self.packing.read.len() * 32).collect_vec();
        let subsequence = &bytes[self.start..self.end];
        if self.same_orientation {
            Vec::from(subsequence)
        } else {
            reverse_complement_ascii_to_ascii(subsequence)
        }
    }
}

impl<Packing> Subsequence<Packing> {
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

#[cfg(debug_assertions)]
fn prefix_str(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    iter_prefix_len(a.iter().copied(), b.iter().copied())
}

#[cfg(debug_assertions)]
fn suffix_str(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    iter_suffix_len(&mut a.iter().copied(), &mut b.iter().copied())
}

fn iter_prefix_len(mut x: impl Iterator<Item = u8>, mut y: impl Iterator<Item = u8>) -> usize {
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        // N is treated like a A
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };
        let yc = if unlikely(yc == b'N') { b'A' } else { yc };

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
        // N is treated like a A
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };
        let yc = if unlikely(yc == b'N') { b'A' } else { yc };

        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}

// #[cfg(test)]
// mod tests {
// use two_bits::encode_2bits;

// use super::*;

// #[test]
// pub fn test_dump() {
//     // TODO
//     // let read = "ACTGAGCTA";
//     // let bytes = read.as_bytes();
//     // let hk = Subsequence::new(bytes, 0, bytes.len(), true);

//     // let dest = &mut [0, 0, 0];
//     // hk.to_canonical().dump_as_2bits(dest);

//     // let expected: Vec<u8> = encode_2bits(bytes.iter().copied(), read.len()).collect();
//     // assert_eq!(dest, expected.as_slice())
// }
// }
