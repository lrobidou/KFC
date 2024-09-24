/// Encodes 32 bases in a `u64`.
#[cfg(target_feature = "bmi2")]
#[target_feature(enable = "bmi2")]
unsafe fn encode_32_bases_bmi(bases: &[u8; 32]) -> u64 {
    use std::arch::x86_64::_pext_u64;

    let bytes_0 = u64::from_be_bytes(*(bases.as_ptr() as *const [u8; 8]));
    let bytes_1 = u64::from_be_bytes(*(bases.as_ptr().add(8) as *const [u8; 8]));
    let bytes_2 = u64::from_be_bytes(*(bases.as_ptr().add(16) as *const [u8; 8]));
    let bytes_3 = u64::from_be_bytes(*(bases.as_ptr().add(24) as *const [u8; 8]));

    unsafe {
        let bytes_0 = _pext_u64(bytes_0, 0x0606060606060606);
        let bytes_1 = _pext_u64(bytes_1, 0x0606060606060606);
        let bytes_2 = _pext_u64(bytes_2, 0x0606060606060606);
        let bytes_3 = _pext_u64(bytes_3, 0x0606060606060606);
        (bytes_0 << 48) | (bytes_1 << 32) | (bytes_2 << 16) | bytes_3
    }
}

/// Encodes up to 32 bases in a `u64`.
/// Likely to be slower than `encode_32_bases`, but available on every architectures.
fn encode_up_to_32_bases_slow(bases: &[u8]) -> u64 {
    assert!(bases.len() <= 32);
    let mut result: u64 = 0;
    for base in bases {
        let base_encoded = ((base >> 1) & 3) as u64;
        result = (result << 2) | base_encoded;
    }
    result << ((32 - bases.len()) * 2)
}

pub fn encode_32_bases_all_arch(bases: &[u8; 32]) -> u64 {
    #[cfg(target_feature = "bmi2")]
    // SAFETY: we are behind a cfg, so it should be fine
    let encoded = unsafe { encode_32_bases_bmi(bases) };
    #[cfg(not(target_feature = "bmi2"))]
    let encoded = encode_up_to_32_bases_slow(bases);
    encoded
}

/// Encodes a sequence of ascii characters in `u64`s.
/// Does not takes into account the reverse complement.
pub struct Encoder<'a> {
    bases: &'a [u8],
    nb_u64_full: usize,
    nb_u64_done: usize,
    is_there_bases_not_full: bool,
}

impl<'a> Encoder<'a> {
    pub fn new(bases: &'a [u8]) -> Self {
        Self {
            bases,
            nb_u64_full: bases.len() / 32,
            nb_u64_done: 0,
            is_there_bases_not_full: bases.len() % 32 != 0,
        }
    }
}

impl<'a> Iterator for Encoder<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        let start = self.nb_u64_done * 32;
        if self.nb_u64_done < self.nb_u64_full {
            let bases = self.bases[start..start + 32].try_into().unwrap();
            let encoded = encode_32_bases_all_arch(bases);
            self.nb_u64_done += 1;
            Some(encoded)
        } else if self.is_there_bases_not_full {
            let rest = encode_up_to_32_bases_slow(&self.bases[start..self.bases.len()]);
            self.is_there_bases_not_full = false;
            Some(rest)
        } else {
            None
        }
    }
}

pub fn encode_32_bases_revcomp(bases: &[u8; 32]) -> u64 {
    let encoded = encode_32_bases_all_arch(bases);
    super::revcomp_32_bases(encoded)
}

pub fn encode_up_to_32_bases_revcomp_slow(bases: &[u8]) -> u64 {
    // this is to be applied before shifting the ASCII encoding
    let mut revcomp_map = [0; 8];
    revcomp_map[0b0000_0000] = 0b0000_0010;
    revcomp_map[0b0000_0100] = 0b0000_0000;
    revcomp_map[0b0000_0010] = 0b0000_0011;
    revcomp_map[0b0000_0110] = 0b0000_0001;

    assert!(bases.len() <= 32);
    let mut result: u64 = 0;
    for base in bases.iter().rev() {
        let base_encoded = (base & 0b0000_0110) as u64;
        let base_encoded_revcomp = revcomp_map[base_encoded as usize];

        result = (result << 2) | base_encoded_revcomp;
    }
    result << ((32 - bases.len()) * 2)
}

pub struct RevCompEncoder<'a> {
    bases: &'a [u8],
    nb_u64_full: usize,
    nb_u64_done: usize,
    nb_bases_not_full: usize,
}

impl<'a> RevCompEncoder<'a> {
    pub fn new(bases: &'a [u8]) -> Self {
        Self {
            bases,
            nb_u64_full: bases.len() / 32,
            nb_u64_done: 0,
            nb_bases_not_full: bases.len() % 32,
        }
    }
}

impl<'a> Iterator for RevCompEncoder<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.nb_u64_done < self.nb_u64_full {
            let end = self.bases.len() - self.nb_u64_done * 32;
            let bases = self.bases[end - 32..end].try_into().unwrap();
            let encoded = encode_32_bases_revcomp(bases);
            self.nb_u64_done += 1;
            Some(encoded)
        } else if self.nb_bases_not_full != 0 {
            let rest = encode_up_to_32_bases_revcomp_slow(&self.bases[0..self.nb_bases_not_full]);
            self.nb_bases_not_full = 0;
            Some(rest)
        } else {
            None
        }
    }
}

// #[target_feature(enable = "bmi2")]
// unsafe fn encode_8_bases(bases: &[u8; 8]) -> u64 {
//     let bytes = u64::from_be_bytes(*bases);
//     unsafe { _pext_u64(bytes, 0x0606060606060606) }
// }

// #[target_feature(enable = "bmi2")]
// unsafe fn encode_16_bases(bases: &[u8]) -> u64 {
//     let bytes_0 = u64::from_be_bytes(*(bases.as_ptr() as *const [u8; 8]));
//     let bytes_1 = u64::from_be_bytes(*(bases.as_ptr().add(8) as *const [u8; 8]));

//     unsafe {
//         let bytes_0 = _pext_u64(bytes_0, 0x0606060606060606);
//         let bytes_1 = _pext_u64(bytes_1, 0x0606060606060606);
//         (bytes_0 << 16) | bytes_1
//     }
// }

/// Encode bases in a vec. Mainly for tests purposes. Prefer using `Encoder`.
#[cfg(test)]
pub fn encode_bases(bases: &[u8]) -> Vec<u64> {
    let nb_u64 = bases.len() / 32;
    let mut v: Vec<u64> = vec![];
    for i in 0..nb_u64 {
        let bases = bases[i * 32..(i + 1) * 32].try_into().unwrap();
        v.push(encode_32_bases_all_arch(bases));
    }
    if bases.len() % 32 != 0 {
        let rest = encode_up_to_32_bases_slow(&bases[nb_u64 * 32..bases.len()]);
        v.push(rest)
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    // >>> bin(ord('A'))
    // '0b01000001'  # 00 <-> 10 0b10111110 =>
    // >>> bin(ord('C'))
    // '0b01000011'  # 01 <-> 11 0b10111100 =>
    // >>> bin(ord('T'))
    // '0b01010100'  # 10 <-> 00 0b10101011 =>
    // >>> bin(ord('G'))
    // '0b01000111'  # 11 <-> 01 0b10111000 =>

    // #[test]
    // fn test_encode_8_bases() {
    //     let read = String::from("TATTTACT");
    //     let mut slice = [0; 8];
    //     slice.copy_from_slice(read.as_bytes());

    //     let encoded = unsafe { encode_8_bases(&slice) };

    //     assert_eq!(
    //         encoded,
    //         // T  A  T  T  T  A  C  T
    //         0b10_00_10_10_10_00_01_10
    //     );
    // }

    // #[test]
    // fn test_encode_16_bases() {
    //     let read = String::from("TATTTACTGTAATGAA");
    //     let read = read.as_bytes();
    //     let encoded = unsafe { encode_16_bases(read) };
    //     assert_eq!(
    //         encoded,
    //         //    TATT     TACT     GTAA     TGAA     GGAC     CTTC     GTCT     CCCC
    //         0b10001010_10000110_11100000_10110000
    //     )
    // }

    #[test]
    fn test_encode_32_bases() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCC");
        let read = read.as_bytes();
        let encoded = encode_32_bases_all_arch(read.try_into().unwrap());
        assert_eq!(
            encoded,
            //    TATT     TACT     GTAA     TGAA     GGAC     CTTC     GTCT     CCCC
            0b10001010_10000110_11100000_10110000_11110001_01101001_11100110_01010101
        )
    }

    #[test]
    fn test_encode_32_bases_revcomp() {
        {
            let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCC");
            let read = read.as_bytes();

            let encoded = encode_32_bases_revcomp(read.try_into().unwrap());
            assert_eq!(
                encoded,
                // G G G G  A G A C  G A A G  G T C C  T T C A  T T A C  A G T A  A A T A
                0b11111111_00110001_11000011_11100101_10100100_10100001_00111000_00001000
            )
        }

        {
            let read = String::from("AAAAAAGGACCTTCGTCTCCCCGGGGAGACGA");
            let read = read.as_bytes();

            let encoded = encode_32_bases_revcomp(read.try_into().unwrap());
            assert_eq!(
                encoded,
                // T C G T  C T C C  C C G G  G G A G  A C G A  A G G T  C C T T  T T T T
                0b10011110_01100101_01011111_11110011_00011100_00111110_01011010_10101010
            )
        }
    }

    #[test]
    fn test_encode_32_bases_revcomp_slow() {
        let read = String::from("TACTGTAATGAAGGACCTTCGTCTCCCC");
        let read = read.as_bytes();

        let encoded = encode_up_to_32_bases_revcomp_slow(read);
        assert_eq!(
            encoded,
            // G G G G  A G A C  G A A G  G T C C  T T C A  T T A C  A G T A
            0b11111111_00110001_11000011_11100101_10100100_10100001_00111000_00000000
        )
    }

    #[test]
    fn test_encode_32_bases_revcomp_slow_unique_letter() {
        let read = String::from("C");
        let read = read.as_bytes();

        let encoded = encode_up_to_32_bases_revcomp_slow(read);
        assert_eq!(
            encoded,
            0b1100000000000000000000000000000000000000000000000000000000000000
        )
    }

    #[test]
    fn test_encode_32_bases_revcomp_slow_3_bases() {
        let read = String::from("CAG");
        let read = read.as_bytes();

        let encoded = encode_up_to_32_bases_revcomp_slow(read);
        assert_eq!(
            encoded,
            0b0110110000000000000000000000000000000000000000000000000000000000
        )
    }

    #[test]
    fn test_encode_32_bases_revcomp_slow_7_bases() {
        let read = String::from("GTTCCAT");
        let read = read.as_bytes();

        let encoded = encode_up_to_32_bases_revcomp_slow(read);
        assert_eq!(
            encoded,
            0b0010111100000100000000000000000000000000000000000000000000000000
        )
    }
}
