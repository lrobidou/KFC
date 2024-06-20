pub fn char_to_uint(c: &u8) -> u64 {
    (*c as u64 >> 1) & 3
}

pub fn uint_to_char(c: &u8) -> u8 {
    match c {
        0 => b'A',
        1 => b'C',
        2 => b'T',
        _ => b'G',
    }
}

// TODO copy
pub fn encode_2bits(bytes: impl Iterator<Item = u8>, len: usize) -> Vec<u64> {
    let add_one = (len % 32) != 0;
    let mut result: Vec<u64> = vec![0; (len / 32) + add_one as usize];
    for (i, ascii_letter) in bytes.into_iter().enumerate() {
        let shift = (31 - i % 32) * 2;
        result[i / 32] += char_to_uint(&ascii_letter) << shift;
    }
    result
}

// TODO copy
pub fn decode_2bits(bytes: impl Iterator<Item = u64>, len: usize) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::with_capacity(len);
    for (i, byte) in bytes.enumerate() {
        if i < len / 32 {
            for j in 0..32 {
                let bits = byte >> (62 - 2 * j) & 0b00000011;
                result.push(uint_to_char(&(bits as u8)));
            }
        } else {
            for j in 0..(len % 32) {
                let bits = byte >> (62 - 2 * j) & 0b00000011;
                result.push(uint_to_char(&(bits as u8)));
            }
        }
    }

    result
}

mod tests {
    use super::*;

    #[test]
    fn test_char_to_uint() {
        let a = b'A';
        let c = b'C';
        let t = b'T';
        let g = b'G';
        assert_eq!(char_to_uint(&a), 0);
        assert_eq!(char_to_uint(&c), 1);
        assert_eq!(char_to_uint(&t), 2);
        assert_eq!(char_to_uint(&g), 3);
    }

    #[test]
    fn test_encode_2bits_multiple_of_4() {
        let seq = "ACGTGCAG";
        // 00011110 11010011
        // = 0001111011010011000000000000000000000000000000000000000000000000
        // 2221119041223786496
        assert_eq!(
            encode_2bits(seq.bytes(), seq.len()),
            Vec::from(&[2221119041223786496])
        );
    }

    #[test]
    fn test_encode_2bits_not_multiple_of_4() {
        let seq = "ACGTGCA";
        // 00011110 11010000 00000000 00000000 00000000 ...  (two last bits unused)
        // = 0001111011010000000000000000000000000000000000000000000000000000
        // 2220274616293654528
        assert_eq!(
            encode_2bits(seq.bytes(), seq.len()),
            Vec::from(&[2220274616293654528])
        );
    }

    #[test]
    fn test_decode_2bits() {
        let seq: [u64; 1] = [0b0001111011010011000000000000000000000000000000000000000000000000];
        assert_eq!(decode_2bits(seq.into_iter(), 8), "ACGTGCAG".as_bytes());
    }

    #[test]
    fn test_back_and_forth() {
        let seq = "ACGTGCAGTAGCATACGACGACATATTAGACAGACATAGACGACTAGACATAGGACATCAGACTATGACGGCAGCATAGCTATTACTCTCTGTGATATACAGTGCATGACTAGTACGACTCAGCATAGCATACGACTAAGCAGCATACGACATCAGACTACGACTACAGCGCGCATTATATTTGCGCTAGCTCCATGAATATAT";
        assert_eq!(
            decode_2bits(encode_2bits(seq.bytes(), seq.len()).into_iter(), seq.len()),
            seq.as_bytes()
        );
    }
}
