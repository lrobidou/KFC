const CHAR_TO_2BIT: [u8; 255] = {
    let mut tab = [0; 255];
    tab[b'A' as usize] = 0;
    tab[b'T' as usize] = 2;
    tab[b'C' as usize] = 1;
    tab[b'G' as usize] = 3;
    tab[b'N' as usize] = 0; // N is A
    tab
};

pub fn char_to_2bit(c: u8) -> u8 {
    unsafe { *CHAR_TO_2BIT.get_unchecked(c as usize) }
}

const TAB_U8_TO_CHAR: [u8; 4] = [b'A', b'C', b'T', b'G'];

pub fn u8_to_char(c: u8) -> u8 {
    unsafe { *TAB_U8_TO_CHAR.get_unchecked(c as usize) }
}

pub fn encode_minimizer(bytes: impl Iterator<Item = u8>) -> u64 {
    const NB_BASES_MAX: usize = 64 / 2;
    let mut nb_base = 0;
    let mut result: u64 = 0;
    for ascii_letter in bytes.into_iter() {
        nb_base += 1;
        result <<= 2;
        result += u64::from(char_to_2bit(ascii_letter));
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

pub fn decode_minimizer(bytes: u64, m: usize) -> String {
    const NB_BASES_MAX: usize = 64 / 2;
    // let mut nb_base = 0;
    // let mut result: u64 = 0;
    let nb_base_not_inserted = NB_BASES_MAX - m;
    let mut bytes = bytes >> (2 * nb_base_not_inserted);

    let mut chars = vec![];
    for _ in 0..m {
        let encoded_base = (bytes & 0b11) as u8;
        bytes >>= 2;
        let char = u8_to_char(encoded_base);
        chars.push(char);
    }
    chars.reverse();
    String::from_utf8(chars).unwrap()
}

// pub fn encode_2bits(bytes: impl Iterator<Item = u8>, len: usize) -> Vec<u8> {
//     let add_one = (len % 4) != 0;
//     let mut result: Vec<u8> = vec![0; (len / 4) + usize::from(add_one)];
//     for (i, ascii_letter) in bytes.into_iter().enumerate() {
//         let shift = (3 - i % 4) * 2;
//         result[i / 4] += char_to_2bit(ascii_letter) << shift;
//     }
//     result
// }

// pub fn encode_2bits<I>(bytes: I, len: usize) -> Encode2Bits<I>
// where
//     I: Iterator<Item = u8>,
// {
//     Encode2Bits {
//         bytes,
//         len,
//         index: 0,
//         buffer: 0,
//         buffer_len: 0,
//     }
// }

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
                self.buffer |= char_to_2bit(ascii_letter) << shift;
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

// pub fn decode_2bits(
//     bytes: &[u8],
//     start: usize,
//     end: usize,
//     nb_base_in_iterator: usize,
// ) -> impl DoubleEndedIterator<Item = u8> + '_ {
//     debug_assert!(end <= nb_base_in_iterator);

//     Decode2Bits {
//         bytes,
//         forward_index: start,
//         backward_index: end,
//     }
// }

pub struct Decode2Bits<'a> {
    bytes: &'a [u8],
    forward_index: usize,
    backward_index: usize,
}

impl<'a> Iterator for Decode2Bits<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.forward_index >= self.backward_index {
            return None;
        }

        let byte_index = self.forward_index / 4;
        let byte = self.bytes[byte_index];
        let shift = 6 - 2 * (self.forward_index % 4);
        let bits = (byte >> shift) & 0b0000_0011;
        let result = u8_to_char(bits);

        self.forward_index += 1;

        Some(result)
    }
}

impl<'a> DoubleEndedIterator for Decode2Bits<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.backward_index <= self.forward_index {
            return None;
        }
        self.backward_index -= 1;
        let byte_index = self.backward_index / 4;
        let byte = self.bytes[byte_index];
        let shift = 6 - 2 * (self.backward_index % 4);
        let bits = (byte >> shift) & 0b0000_0011;
        let result = u8_to_char(bits);

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_char_to_u64() {
        let a = b'A';
        let c = b'C';
        let t = b'T';
        let g = b'G';
        assert_eq!(char_to_2bit(a), 0);
        assert_eq!(char_to_2bit(c), 1);
        assert_eq!(char_to_2bit(t), 2);
        assert_eq!(char_to_2bit(g), 3);
    }

    #[test]
    fn encode_decode() {
        let minimizer = String::from("ACTGCAGCTA");
        let m = minimizer.len();

        let encoded = encode_minimizer(minimizer.as_bytes().iter().copied());
        let decoded = decode_minimizer(encoded, m);

        assert_eq!(decoded, minimizer);
    }
}
