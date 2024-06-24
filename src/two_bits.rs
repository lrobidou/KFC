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

// get unchecked
pub fn u8_to_char(c: u8) -> u8 {
    match c {
        0 => b'A',
        1 => b'C',
        2 => b'T',
        _ => b'G',
    }
}

// TODO copy
pub fn encode_2bits_u64(bytes: impl Iterator<Item = u8>, len: usize) -> Vec<u64> {
    let add_one = (len % 32) != 0;
    let mut result: Vec<u64> = vec![0; (len / 32) + usize::from(add_one)];
    for (i, ascii_letter) in bytes.into_iter().enumerate() {
        let shift = (31 - i % 32) * 2;
        result[i / 32] += char_to_u64(ascii_letter) << shift;
    }
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

// // TODO copy
// pub fn decode_2bits(bytes: impl Iterator<Item = u64>, len: usize) -> Vec<u8> {
//     let mut result: Vec<u8> = Vec::with_capacity(len);
//     for (i, byte) in bytes.enumerate() {
//         if i < len / 32 {
//             for j in 0..32 {
//                 let bits = byte >> (62 - 2 * j) & 0b0000_0011;
//                 result.push(u8_to_char(bits as u8));
//             }
//         } else {
//             for j in 0..(len % 32) {
//                 let bits = byte >> (62 - 2 * j) & 0b0000_0011;
//                 result.push(u8_to_char(bits as u8));
//             }
//         }
//     }

//     result
// }

// // TODO copy
// pub fn decode_2bits(bytes: impl Iterator<Item = u8>, len: usize) -> Vec<u8> {
//     let mut result: Vec<u8> = Vec::with_capacity(len);
//     for (i, byte) in bytes.enumerate() {
//         if i < len / 4 {
//             for j in 0..4 {
//                 let bits = byte >> (6 - 2 * j) & 0b0000_0011;
//                 result.push(u8_to_char(bits));
//             }
//         } else {
//             for j in 0..(len % 4) {
//                 let bits = byte >> (6 - 2 * j) & 0b0000_0011;
//                 result.push(u8_to_char(bits));
//             }
//         }
//     }

//     result
// }

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

    let nb_base_in_last_byte = nb_base_in_iterator % 4;
    let nb_full_bases = nb_base_in_iterator - nb_base_in_last_byte;
    // TODO relire
    let full_bytes_to_skip_at_the_end = if end > nb_full_bases {
        assert!(end < nb_full_bases + 4);
        0
    } else if end == nb_full_bases {
        if nb_base_in_last_byte != 0 {
            1
        } else {
            0
        }
    } else {
        let to_keep = (end + 3) & !3; // equivalent to `let to_keep = 4 * ((end / 4) + (end % 4 != 0) as usize);`
        let skip = nb_full_bases - to_keep;
        skip / 4 + 1
    };

    for _ in 0..full_bytes_to_skip_at_the_end {
        bytes.next_back();
    }

    Decode2Bits {
        bytes,
        forward_index: start,
        byte: 0,
        byte_ready: false,
        end_byte: None,
        backward_index: end,
        end_byte_ready: false,
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
            // print!("no bytes ready");
            self.byte = self.bytes.next()?;
            self.byte_ready = true;
        } else {
            // print!(".");
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
        // println!("backward index: {}", self.backward_index);
        if self.backward_index <= self.forward_index {
            return None;
        }

        if !self.end_byte_ready {
            // println!("pulling a byte");
            self.end_byte = self.bytes.next_back();
            self.end_byte_ready = true;
        } else {
            // print!(".");
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
