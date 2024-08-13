pub fn char_to_u64(c: u8) -> u64 {
    (u64::from(c) >> 1) & 3
}

// // get unchecked
// pub fn u64_to_char(c: u8) -> u8 {
//     match c {
//         0 => b'A',
//         1 => b'C',
//         2 => b'T',
//         _ => b'G',
//     }
// }

pub fn char_to_u8(c: u8) -> u8 {
    (c >> 1) & 3
}

const TAB_U8_TO_CHAR: [u8; 4] = [b'A', b'C', b'T', b'G'];

pub fn u8_to_char(c: u8) -> u8 {
    unsafe { *TAB_U8_TO_CHAR.get_unchecked(c as usize) }
}

pub fn encode_minimizer(bytes: impl Iterator<Item = u8>) -> u64 {
    const NB_BASES_MAX: usize = 64 / 2;
    let mut nb_base = 0;
    let mut result = 0;
    for ascii_letter in bytes.into_iter() {
        nb_base += 1;
        result <<= 2;
        result += char_to_u64(ascii_letter);
    }
    assert!(
        nb_base <= NB_BASES_MAX,
        "minimizer should consists of {NB_BASES_MAX} bases max"
    );
    // TODO possible source of bug if kff need some other order
    let nb_base_not_inserted = NB_BASES_MAX - nb_base;
    let result = result << (2 * nb_base_not_inserted);
    #[allow(clippy::let_and_return)]
    result
}

// pub fn encode_2bits(bytes: impl Iterator<Item = u8>, len: usize) -> Vec<u8> {
//     let add_one = (len % 4) != 0;
//     let mut result: Vec<u8> = vec![0; (len / 4) + usize::from(add_one)];
//     for (i, ascii_letter) in bytes.into_iter().enumerate() {
//         let shift = (3 - i % 4) * 2;
//         result[i / 4] += char_to_u8(ascii_letter) << shift;
//     }
//     result
// }

pub fn encode_2bits<I>(bytes: I, len: usize) -> Encode2Bits<I>
where
    I: Iterator<Item = u8>,
{
    Encode2Bits {
        bytes,
        len,
        index: 0,
        buffer: 0,
        buffer_len: 0,
    }
}

pub struct Encode2Bits<I>
where
    I: Iterator<Item = u8>,
{
    bytes: I,
    len: usize,
    index: usize,
    buffer: u8,
    buffer_len: usize,
}

impl<I> Iterator for Encode2Bits<I>
where
    I: Iterator<Item = u8>,
{
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.len {
            return None;
        }

        let mut result = 0;
        while self.buffer_len < 8 && self.index < self.len {
            if let Some(ascii_letter) = self.bytes.next() {
                let shift = (3 - self.index % 4) * 2;
                self.buffer |= char_to_u8(ascii_letter) << shift;
                self.buffer_len += 2;
                self.index += 1;
            } else {
                break;
            }
        }

        if self.buffer_len >= 8 {
            result = self.buffer;
            self.buffer = 0;
            self.buffer_len -= 8;
        } else if self.buffer_len > 0 {
            result = self.buffer;
            self.buffer = 0;
            self.buffer_len = 0;
        }

        Some(result)
    }
}

pub fn decode_2bits<I>(
    mut bytes: I,
    start: usize,
    end: usize,
    nb_base_in_iterator: usize,
) -> impl DoubleEndedIterator<Item = u8>
where
    I: DoubleEndedIterator<Item = u8>,
{
    debug_assert!(end <= nb_base_in_iterator);
    // Calculate how many full bytes to skip at the start
    let full_bytes_to_skip = start / 4;
    for _ in 0..full_bytes_to_skip {
        bytes.next();
    }

    // TODO relire
    let nb_base_in_last_byte = nb_base_in_iterator % 4;
    let nb_full_bases: usize = nb_base_in_iterator - nb_base_in_last_byte;

    Decode2Bits {
        bytes,
        forward_index: start,
        byte: 0,
        byte_ready: false,
        end_byte: None,
        backward_index: end,
        end_byte_ready: false,
        iterating_back: false,
        nb_base_in_last_byte,
        nb_full_bases,
    }
}

pub struct Decode2Bits<I>
where
    I: DoubleEndedIterator<Item = u8>,
{
    bytes: I,
    forward_index: usize,
    byte: u8,
    byte_ready: bool,
    end_byte: Option<u8>,
    backward_index: usize,
    end_byte_ready: bool,
    iterating_back: bool,
    nb_base_in_last_byte: usize,
    nb_full_bases: usize,
}

impl<I> Iterator for Decode2Bits<I>
where
    I: DoubleEndedIterator<Item = u8>,
{
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.forward_index >= self.backward_index {
            return None;
        }

        if !self.byte_ready {
            self.byte = self.bytes.next()?;
            self.byte_ready = true;
        }
        let shift = 6 - 2 * (self.forward_index % 4);
        let bits = (self.byte >> shift) & 0b0000_0011;
        let result = u8_to_char(bits);

        self.forward_index += 1;
        if self.forward_index % 4 == 0 {
            self.byte_ready = false;
        }

        Some(result)
    }
}

impl<I> DoubleEndedIterator for Decode2Bits<I>
where
    I: DoubleEndedIterator<Item = u8>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        if !self.iterating_back {
            let end = self.backward_index;

            #[allow(clippy::comparison_chain)] // I find it clearer than the alternative
            let full_bytes_to_skip_at_the_end = if end > self.nb_full_bases {
                assert!(end < self.nb_full_bases + 4);
                0
            } else if end == self.nb_full_bases {
                if self.nb_base_in_last_byte != 0 {
                    1
                } else {
                    0
                }
            } else {
                let to_keep = (end + 3) & !3; // equivalent to `let to_keep = 4 * ((end / 4) + (end % 4 != 0) as usize);`

                let skip = self.nb_full_bases - to_keep;
                skip / 4 + 1
            };

            for _ in 0..full_bytes_to_skip_at_the_end {
                self.bytes.next_back();
            }
            self.iterating_back = true;
        }

        if self.backward_index <= self.forward_index {
            return None;
        }

        if !self.end_byte_ready {
            self.end_byte = self.bytes.next_back();
            self.end_byte_ready = true;
        }

        let shift = 6 - 2 * ((self.backward_index - 1) % 4);
        let bits = (self.end_byte? >> shift) & 0b0000_0011;
        let result = u8_to_char(bits);

        self.backward_index -= 1;
        if self.backward_index % 4 == 0 {
            self.end_byte_ready = false;
        }

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_char_to_u64() {
        let a = b'A';
        let c = b'C';
        let t = b'T';
        let g = b'G';
        assert_eq!(char_to_u8(a), 0);
        assert_eq!(char_to_u8(c), 1);
        assert_eq!(char_to_u8(t), 2);
        assert_eq!(char_to_u8(g), 3);
    }

    #[test]
    fn test_encode_2bits_multiple_of_4() {
        let seq = "ACGTGCAG";
        // 00011110 11010011
        // = 0001111011010011000000000000000000000000000000000000000000000000
        // 2221119041223786496
        assert_eq!(
            encode_2bits(seq.bytes(), seq.len()).collect::<Vec<_>>(),
            Vec::from(&[0b00011110, 0b11010011])
        );
    }

    #[test]
    fn test_encode_2bits_not_multiple_of_4() {
        let seq = "ACGTGCA";
        // 00011110 11010000 00000000 00000000 00000000 ...  (two last bits unused)
        // = 0001111011010000000000000000000000000000000000000000000000000000
        // 2220274616293654528
        assert_eq!(
            encode_2bits(seq.bytes(), seq.len()).collect::<Vec<_>>(),
            Vec::from(&[0b00011110, 0b11010000])
        );
    }

    #[test]
    fn test_decode_2bits() {
        // let seq: [u64; 1] = [0b0001111011010011000000000000000000000000000000000000000000000000];
        let seq: [u8; 2] = [0b00011110, 0b11010011];
        assert_eq!(
            decode_2bits(seq.into_iter(), 0, 8, 8).collect_vec(),
            "ACGTGCAG".as_bytes()
        );
    }

    #[test]
    fn test_decode_2bits_beginning_not_aligned() {
        // let seq: [u64; 1] = [0b0001111011010011000000000000000000000000000000000000000000000000];
        let seq: [u8; 2] = [0b00011110, 0b11010011];
        assert_eq!(
            decode_2bits(seq.into_iter(), 3, 8, 8).collect_vec(),
            "TGCAG".as_bytes()
        );
    }

    #[test]
    fn test_encode_decode() {
        let seq = "ACGTGCAGTAGCATACGACGACATATTAGACAGACATAGACGACTAGACATAGGACATCAGACTATGACGGCAGCATAGCTATTACTCTCTGTGATATACAGTGCATGACTAGTACGACTCAGCATAGCATACGACTAAGCAGCATACGACATCAGACTACGACTACAGCGCGCATTATATTTGCGCTAGCTCCATGAATATAT";
        assert_eq!(
            decode_2bits(
                encode_2bits(seq.bytes(), seq.len())
                    .collect_vec()
                    .into_iter(),
                0,
                seq.len(),
                seq.len()
            )
            .collect_vec(),
            seq.as_bytes()
        );
    }

    #[test]
    fn test_decode_reverse() {
        let seq = "ACGTTGCAG";
        let encode = encode_2bits(seq.bytes(), seq.len()).collect_vec();

        let decode =
            decode_2bits(encode.clone().into_iter(), 0, seq.len(), seq.len()).collect_vec();
        let decode_rev = decode_2bits(encode.into_iter(), 0, seq.len(), seq.len())
            .rev()
            .collect_vec()
            .iter()
            .rev()
            .copied()
            .collect_vec();

        assert_eq!(decode, decode_rev);
    }
}
