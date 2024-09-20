mod bmi;
pub use bmi::{Decoder, Encoder};

use std::arch::x86_64::{
    __m256i, _mm256_and_si256, _mm256_loadu_si256, _mm256_or_si256, _mm256_permute2f128_si256,
    _mm256_set1_epi8, _mm256_setr_epi8, _mm256_shuffle_epi8, _mm256_slli_epi16, _mm256_srli_epi16,
    _mm256_xor_si256, _pext_u64,
};

// >>> bin(ord('A'))
// '0b01000001'  # 00 <-> 10 0b10111110 =>
// >>> bin(ord('C'))
// '0b01000011'  # 01 <-> 11 0b10111100 =>
// >>> bin(ord('T'))
// '0b01010100'  # 10 <-> 00 0b10101011 =>
// >>> bin(ord('G'))
// '0b01000111'  # 11 <-> 01 0b10111000 =>

// adapted from an idea found here https://fmilicchio.github.io/ymm0.html
// changes:
//  - use AVX2
//  - different computation of the encoding
/// Takes a slice ASCII chars ({A, C, T, G}) and computes the 2 bits encoding of the DNA string.
/// Result is undefined for non-ACTG chars.
/// SAFETY: there must be at least 128 bases in the slice.
unsafe fn encode_128_bases(bases: &[u8]) -> __m256i {
    assert!(is_x86_feature_detected!("avx2"));

    let ptr = bases.as_ptr() as *const __m256i;

    unsafe {
        // latency, throughput = 7, 0.56
        let register_0 = _mm256_loadu_si256(ptr);
        let register_1 = _mm256_loadu_si256(ptr.add(1));
        let register_2 = _mm256_loadu_si256(ptr.add(2));
        let register_3 = _mm256_loadu_si256(ptr.add(3));
        // total: 8.5 cycles (+8.5 cycles)

        // mask to get the two least significant bits of each byte
        let mask = _mm256_set1_epi8(0b0000_0011);

        // shift by one bit (as if it's a vector of *16* bits)
        // latency, throughput = 1, 0.5
        let register_0 = _mm256_srli_epi16(register_0, 1);
        let register_1 = _mm256_srli_epi16(register_1, 1);
        let register_2 = _mm256_srli_epi16(register_2, 1);
        let register_3 = _mm256_srli_epi16(register_3, 1);
        // total: 11 cycles (+2.5 cycles)

        // get the two least significant bits of each byte
        // latency, throughput = 1, 0.333333333
        let register_0 = _mm256_and_si256(register_0, mask);
        let register_1 = _mm256_and_si256(register_1, mask);
        let register_2 = _mm256_and_si256(register_2, mask);
        let register_3 = _mm256_and_si256(register_3, mask);
        // total: 13 cycles (+2 cycles)

        // shift the data so that two consecutive indexes do not overlap
        // latency, throughput = 1, 0.5
        // let register_0 = _mm256_slli_epi16::<0>(register_0);
        let register_1 = _mm256_slli_epi16(register_1, 2);
        let register_2 = _mm256_slli_epi16(register_2, 4);
        let register_3 = _mm256_slli_epi16(register_3, 6);
        // total: 15 cycles (+2 cycles)

        // AND everything
        // latency, throughput = 1, 0.333333333
        let register_0 = _mm256_or_si256(register_0, register_1);
        let register_2 = _mm256_or_si256(register_2, register_3);

        _mm256_or_si256(register_0, register_2)
        // total: 17.3 cycles (+2.3 cycles)
    }
}

// /// Computes the ascii ogf the string in the register.
// /// The lookup table should be empty, apart from the 4 first bytes.
// unsafe fn to_ascii_128_bases(register: __m256i, lookup_table: &[u8; 256]) {
//     // Load the first 32 bytes, but only the first 4 will be non-zero in this case
//     let lookup_chunk = _mm256_loadu_si256(lookup_table.as_ptr() as *const __m256i);

//     // masks to get the two bits of each byte
//     let mask_0 = _mm256_set1_epi8(0b0000_0011);
//     let mask_1 = _mm256_set1_epi8(0b0000_1100);
//     let mask_2 = _mm256_set1_epi8(0b0011_0000);
//     let mask_3 = _mm256_set1_epi8(0b1100_0000u8 as i8);

//     // get the two least significant bits of each byte
//     // latency, throughput = 1, 0.333333333
//     let register_0 = _mm256_and_si256(register, mask_0);
//     let register_1 = _mm256_and_si256(register, mask_1);
//     let register_2 = _mm256_and_si256(register, mask_2);
//     let register_3 = _mm256_and_si256(register, mask_3);

//     // shift the data so that indexes are aligned
//     // latency, throughput = 1, 0.5
//     // let register_0 = _mm256_slli_epi16::<0>(register_0);
//     let register_1 = _mm256_srli_epi16(register_1, 2);
//     let register_2 = _mm256_srli_epi16(register_2, 4);
//     let register_3 = _mm256_srli_epi16(register_3, 6);

//     // Lookup using shuffle for only the first chunk (which includes the first 4 non-zero values)
//     let register_0 = _mm256_shuffle_epi8(lookup_chunk, register_0);
//     let register_1 = _mm256_shuffle_epi8(lookup_chunk, register_1);
//     let register_2 = _mm256_shuffle_epi8(lookup_chunk, register_2);
//     let register_3 = _mm256_shuffle_epi8(lookup_chunk, register_3);

//     // Any value greater than or equal to 4 should become zero
//     let is_gte_4 = _mm256_cmpgt_epi8(input, _mm256_set1_epi8(3)); // Compare input > 3
//     let zeroed_result = _mm256_blendv_epi8(result, _mm256_setzero_si256(), is_gte_4);

//     zeroed_result
// }

/// Takes a slice ASCII chars ({A, C, T, G}) and computes the 2 bits encoding of the DNA string.
/// Result is undefined for non-ACTG chars.
/// SAFETY: there must be at least 128 bases in the slice.
unsafe fn encode_128_bases_revcomp(bases: &[u8]) -> __m256i {
    assert!(is_x86_feature_detected!("avx2"));

    let ptr = bases.as_ptr() as *const __m256i;
    unsafe {
        // mask to get the two least significant bits of each byte
        let mask = _mm256_set1_epi8(0b0000_0011);

        // mask to reverse the bytes of a register
        let shuffle_mask = _mm256_setr_epi8(
            31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
            9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        );

        // mask to get the odd bits of a register
        let odd_mask = _mm256_set1_epi8(0b1010_1010u8 as i8);

        // latency, throughput = 7, 0.56
        let register_3 = _mm256_loadu_si256(ptr);
        let register_2 = _mm256_loadu_si256(ptr.add(1));
        let register_1 = _mm256_loadu_si256(ptr.add(2));
        let register_0 = _mm256_loadu_si256(ptr.add(3));
        // total: 8.5 cycles (+8.5 cycles)

        // reverse the bytes of the registers
        let register_0 = _mm256_shuffle_epi8(register_0, shuffle_mask);
        let register_1 = _mm256_shuffle_epi8(register_1, shuffle_mask);
        let register_2 = _mm256_shuffle_epi8(register_2, shuffle_mask);
        let register_3 = _mm256_shuffle_epi8(register_3, shuffle_mask);

        let register_0 = _mm256_permute2f128_si256(register_0, register_0, 1);
        let register_1 = _mm256_permute2f128_si256(register_1, register_1, 1);
        let register_2 = _mm256_permute2f128_si256(register_2, register_2, 1);
        let register_3 = _mm256_permute2f128_si256(register_3, register_3, 1);

        // shift by one bit (as if it's a vector of *16* bits)
        // latency, throughput = 1, 0.5
        let register_0 = _mm256_srli_epi16(register_0, 1);
        let register_1 = _mm256_srli_epi16(register_1, 1);
        let register_2 = _mm256_srli_epi16(register_2, 1);
        let register_3 = _mm256_srli_epi16(register_3, 1);
        // total: 11 cycles (+2.5 cycles)

        // get the two least significant bits of each byte
        // latency, throughput = 1, 0.333333333
        let register_0 = _mm256_and_si256(register_0, mask);
        let register_1 = _mm256_and_si256(register_1, mask);
        let register_2 = _mm256_and_si256(register_2, mask);
        let register_3 = _mm256_and_si256(register_3, mask);
        // total: 13 cycles (+2 cycles)

        // shift the data so that two consecutive indexes do not overlap
        // let register_0 = _mm256_slli_epi16::<0>(register_0);
        // latency, throughput = 1, 0.5
        let register_1 = _mm256_slli_epi16::<2>(register_1);
        let register_2 = _mm256_slli_epi16::<4>(register_2);
        let register_3 = _mm256_slli_epi16::<6>(register_3);
        // total: 15 cycles (+2 cycles)

        // OR everything
        // latency, throughput = 1, 0.333333333
        let register_0 = _mm256_or_si256(register_0, register_1);
        let register_2 = _mm256_or_si256(register_2, register_3);

        let register = _mm256_or_si256(register_0, register_2);
        // total: 17.3 cycles (+2.3 cycles)

        // XOR the input with the odd mask to flip the odd bits
        _mm256_xor_si256(register, odd_mask)
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// with small functions as intermediate steps
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Takes a pointer to a slice of ASCII chars ({A, C, T, G}) and computes the 2 bits encoding of the DNA string. Then shift it.
/// Result is undefined for non-ACTG chars.
/// SAFETY: there must be at least 32*`idx` bases in the slice.
unsafe fn encode_register<const SHIFT: i32>(register: __m256i) -> __m256i {
    unsafe {
        // mask to get the two least significant bits of each byte
        let mask = _mm256_set1_epi8(0b0000_0011);

        // shift by one bit (as if it's a vector of *16* bits)
        let register: __m256i = _mm256_srli_epi16(register, 1);
        // get the two least significant bits of each byte
        let register = _mm256_and_si256(register, mask);
        // shift the data so that two consecutive indexes do not overlap
        let register = _mm256_slli_epi16::<SHIFT>(register);

        #[allow(clippy::let_and_return)]
        register
    }
}

/// Takes a slice ASCII chars ({A, C, T, G}) and computes the 2 bits encoding of the DNA string.
/// Result is undefined for non-ACTG chars.
/// SAFETY: there must be at least 128 bases in the slice.
unsafe fn encode_128_bases_slow(bases: &[u8]) -> __m256i {
    assert!(is_x86_feature_detected!("avx2"));

    let ptr = bases.as_ptr() as *const __m256i;

    unsafe {
        _mm256_or_si256(
            _mm256_or_si256(
                encode_register::<0>(_mm256_loadu_si256(ptr.add(0))),
                encode_register::<2>(_mm256_loadu_si256(ptr.add(1))),
            ),
            _mm256_or_si256(
                encode_register::<4>(_mm256_loadu_si256(ptr.add(2))),
                encode_register::<6>(_mm256_loadu_si256(ptr.add(3))),
            ),
        )
    }
}

unsafe fn reverse_order_bytes(register: __m256i) -> __m256i {
    // safe as I'm passing valid bytes
    unsafe {
        let shuffle_mask = _mm256_setr_epi8(
            31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
            9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
        );
        // TODO prove safe
        // once proven save, remove the unsafe from the function definition
        let register = _mm256_shuffle_epi8(register, shuffle_mask);
        _mm256_permute2f128_si256(register, register, 1)
    }
}

/// Computes the complement of the bases in `register`.
fn complement(register: __m256i) -> __m256i {
    // strategy:
    // In register, we have the folloowing two_bits encoding
    // A <=> 00
    // C <=> 01
    // T <=> 10
    // G <=> 11
    // the reverse complements are therefore 00 <=> 10 and 01 <=> 11
    // notice we can keep the last bit and complement the other one
    unsafe {
        // mask to get the two least significant bits of each byte
        let odd_mask: __m256i = _mm256_set1_epi8(0b1010_1010u8 as i8);
        // XOR the input with the odd mask to flip the odd bits
        _mm256_xor_si256(register, odd_mask)
    }
}

/// Takes a slice ASCII chars ({A, C, T, G}) and computes the 2 bits encoding of the DNA string.
/// Result is undefined for non-ACTG chars.
/// SAFETY: there must be at least 128 bases in the slice.
fn encode_128_bases_revcomp_slow(bases: &[u8]) -> __m256i {
    assert!(is_x86_feature_detected!("avx2"));

    let ptr = bases.as_ptr() as *const __m256i;

    unsafe {
        complement(_mm256_and_si256(
            _mm256_and_si256(
                encode_register::<0>(reverse_order_bytes(_mm256_loadu_si256(ptr.add(0)))),
                encode_register::<2>(reverse_order_bytes(_mm256_loadu_si256(ptr.add(1)))),
            ),
            _mm256_and_si256(
                encode_register::<4>(reverse_order_bytes(_mm256_loadu_si256(ptr.add(2)))),
                encode_register::<6>(reverse_order_bytes(_mm256_loadu_si256(ptr.add(3)))),
            ),
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::arch::x86_64::{_mm256_cmpeq_epi8, _mm256_movemask_epi8};

    use super::*;
    // >>> bin(ord('A'))
    // '0b01000001'  # 00 <-> 10 0b10111110 =>
    // >>> bin(ord('C'))
    // '0b01000011'  # 01 <-> 11 0b10111100 =>
    // >>> bin(ord('T'))
    // '0b01010100'  # 10 <-> 00 0b10101011 =>
    // >>> bin(ord('G'))
    // '0b01000111'  # 11 <-> 01 0b10111000 =>
    ////////////////////////////////////////////////////////////////////////////////

    fn are_equal(a: __m256i, b: __m256i) -> bool {
        unsafe {
            // Compare a and b for equality, byte by byte
            let cmp = _mm256_cmpeq_epi8(a, b);

            // Create a mask: each bit is set if the corresponding byte is 0xFF
            let mask = _mm256_movemask_epi8(cmp);

            // If all 32 bytes are equal, mask will be -1 (all bits set)
            mask == -1
        }
    }

    fn reinterpret_as_u8_slice(register: __m256i) -> [u8; 32] {
        unsafe {
            // Cast the __m256i as a pointer to a 32-element array of u8
            *(&register as *const __m256i as *const [u8; 32])
        }
    }

    fn decode_register(register: __m256i) -> [u8; 128] {
        let slice = reinterpret_as_u8_slice(register);
        let mut decode = [0; 128];
        for i in 0..128 {
            let byte = slice[i];
            for j in 0..4 {
                let shift = 3 - j;
                let bits = (byte >> shift) & 0b0000_0011;
                decode[i] = match bits {
                    0 => b'A',
                    1 => b'C',
                    2 => b'T',
                    3 => b'G',
                    _ => unreachable!(),
                };
            }
        }
        decode
    }

    #[test]
    fn test_reverse_order_bytes() {
        unsafe {
            let register = _mm256_setr_epi8(
                0, 4, 8, 12, 16, 20, 24, 28, //
                1, 5, 9, 13, 17, 21, 25, 29, //
                2, 6, 10, 14, 18, 22, 26, 30, //
                3, 7, 11, 15, 19, 23, 27, 31, //
            );
            let register_reverse = _mm256_setr_epi8(
                31, 27, 23, 19, 15, 11, 7, 3, //
                30, 26, 22, 18, 14, 10, 6, 2, //
                29, 25, 21, 17, 13, 9, 5, 1, //
                28, 24, 20, 16, 12, 8, 4, 0,
            );
            assert!(are_equal(reverse_order_bytes(register), register_reverse));
        }
    }

    #[test]
    fn test_complement() {
        // A <=> 00
        // C <=> 01
        // T <=> 10
        // G <=> 11
        let mut bases: [u8; 32] = [0b0000_0000; 32];
        bases[0] = 0b0001_1011;
        bases[2] = 0b1111_1111;

        let mut the_complement_i_need: [u8; 32] = [0b1010_1010; 32];
        the_complement_i_need[0] = 0b1011_0001;
        the_complement_i_need[2] = 0b0101_0101;

        let ptr = bases.as_ptr() as *const __m256i;

        let the_complement_i_deserve = unsafe {
            let register = _mm256_loadu_si256(ptr);
            let the_complement_i_deserve = complement(register);
            reinterpret_as_u8_slice(the_complement_i_deserve)
        };

        // you get what you deserve... hopefully
        assert_eq!(the_complement_i_need, the_complement_i_deserve)
    }

    #[test]
    fn test_encode_register() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCC");
        let read = read.as_bytes();
        let mut slice = [0; 32];
        slice.copy_from_slice(read);
        let ptr = slice.as_ptr() as *const __m256i;
        unsafe {
            let register = _mm256_loadu_si256(ptr);

            // shift 0
            let result_shift_0 = encode_register::<0>(register);
            let result_shift_0 = reinterpret_as_u8_slice(result_shift_0);
            assert_eq!(
                result_shift_0,
                [
                    0b00000010, 0b00000000, 0b00000010, 0b00000010, 0b00000010, 0b00000000,
                    0b00000001, 0b00000010, 0b00000011, 0b00000010, 0b00000000, 0b00000000,
                    0b00000010, 0b00000011, 0b00000000, 0b00000000, 0b00000011, 0b00000011,
                    0b00000000, 0b00000001, 0b00000001, 0b00000010, 0b00000010, 0b00000001,
                    0b00000011, 0b00000010, 0b00000001, 0b00000010, 0b00000001, 0b00000001,
                    0b00000001, 0b00000001,
                ]
            );

            // shift 2
            let result_shift_2 = encode_register::<2>(register);
            let result_shift_2 = reinterpret_as_u8_slice(result_shift_2);
            assert_eq!(
                result_shift_2,
                [
                    0b00001000, 0b00000000, 0b00001000, 0b00001000, 0b00001000, 0b00000000,
                    0b00000100, 0b00001000, 0b00001100, 0b00001000, 0b00000000, 0b00000000,
                    0b00001000, 0b00001100, 0b00000000, 0b00000000, 0b00001100, 0b00001100,
                    0b00000000, 0b00000100, 0b00000100, 0b00001000, 0b00001000, 0b00000100,
                    0b00001100, 0b00001000, 0b00000100, 0b00001000, 0b00000100, 0b00000100,
                    0b00000100, 0b00000100,
                ]
            );

            // shift 4
            let result_shift_4 = encode_register::<4>(register);
            let result_shift_4 = reinterpret_as_u8_slice(result_shift_4);
            assert_eq!(
                result_shift_4,
                [
                    0b00100000, 0b00000000, 0b00100000, 0b00100000, 0b00100000, 0b00000000,
                    0b00010000, 0b00100000, 0b00110000, 0b00100000, 0b00000000, 0b00000000,
                    0b00100000, 0b00110000, 0b00000000, 0b00000000, 0b00110000, 0b00110000,
                    0b00000000, 0b00010000, 0b00010000, 0b00100000, 0b00100000, 0b00010000,
                    0b00110000, 0b00100000, 0b00010000, 0b00100000, 0b00010000, 0b00010000,
                    0b00010000, 0b00010000,
                ]
            );

            // shift 6
            let result_shift_6 = encode_register::<6>(register);
            let result_shift_6 = reinterpret_as_u8_slice(result_shift_6);
            assert_eq!(
                result_shift_6,
                [
                    0b10000000, 0b00000000, 0b10000000, 0b10000000, 0b10000000, 0b00000000,
                    0b01000000, 0b10000000, 0b11000000, 0b10000000, 0b00000000, 0b00000000,
                    0b10000000, 0b11000000, 0b00000000, 0b00000000, 0b11000000, 0b11000000,
                    0b00000000, 0b01000000, 0b01000000, 0b10000000, 0b10000000, 0b01000000,
                    0b11000000, 0b10000000, 0b01000000, 0b10000000, 0b01000000, 0b01000000,
                    0b01000000, 0b01000000,
                ]
            );
        }
    }

    #[test]
    fn test_encode_string_slow() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCTGGGCGAGGGGCCAAGGTTGCGCTTTCAGAGCGTTAATATATCCAGTGTAATCAATGATATAATATGGACAGATCATGAAGGTAACTGTTCTTTAG");
        let read = read.as_bytes();
        let mut slice = [0; 128];
        slice.copy_from_slice(read);

        unsafe {
            let register = encode_128_bases_slow(&slice);

            let register = reinterpret_as_u8_slice(register);
            assert_eq!(
                register,
                [
                    // TATTTACTGTAATGAAGGACCTTCGTCTCCCC
                    // TGGGCGAGGGGCCAAGGTTGCGCTTTCAGAGC
                    // GTTAATATATCCAGTGTAATCAATGATATAAT
                    // ATGGACAGATCATGAAGGTAACTGTTCTTTAG

                    // look at the first column: TTGA
                    // then look at the first byte: AGTT
                    // it's storesd in reverse
                    0b00111010, 0b10101100, 0b11101110, 0b11001110, 0b00000110, 0b01101100,
                    0b00000001, 0b11101110, 0b00001111, 0b10101110, 0b01011100, 0b00010100,
                    0b10000110, 0b11110011, 0b00100000, 0b00111100, 0b11101111, 0b11001011,
                    0b10001000, 0b00101101, 0b00010101, 0b01001110, 0b10000110, 0b11101001,
                    0b10111011, 0b10001010, 0b01100101, 0b10000010, 0b10101101, 0b10000001,
                    0b00001101, 0b11100101,
                ]
            );
        }
    }

    #[test]
    fn test_encode_string_fast() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCTGGGCGAGGGGCCAAGGTTGCGCTTTCAGAGCGTTAATATATCCAGTGTAATCAATGATATAATATGGACAGATCATGAAGGTAACTGTTCTTTAG");
        let read = read.as_bytes();
        let mut slice = [0; 128];
        slice.copy_from_slice(read);

        unsafe {
            let register = encode_128_bases(&slice);

            let register = reinterpret_as_u8_slice(register);
            assert_eq!(
                register,
                [
                    // divide the string in 4:
                    // TATTTACTGTAATGAAGGACCTTCGTCTCCCC
                    // TGGGCGAGGGGCCAAGGTTGCGCTTTCAGAGC
                    // GTTAATATATCCAGTGTAATCAATGATATAAT
                    // ATGGACAGATCATGAAGGTAACTGTTCTTTAG
                    //
                    // look at the first column: TTGA
                    // then look at the first byte: 0b00111010 <=> AGTT
                    //
                    // each byte represent a column of 4 letters
                    0b00111010, 0b10101100, 0b11101110, 0b11001110, 0b00000110, 0b01101100,
                    0b00000001, 0b11101110, 0b00001111, 0b10101110, 0b01011100, 0b00010100,
                    0b10000110, 0b11110011, 0b00100000, 0b00111100, 0b11101111, 0b11001011,
                    0b10001000, 0b00101101, 0b00010101, 0b01001110, 0b10000110, 0b11101001,
                    0b10111011, 0b10001010, 0b01100101, 0b10000010, 0b10101101, 0b10000001,
                    0b00001101, 0b11100101,
                ]
            );
        }
    }

    #[test]
    fn test_encode_string_revcomp_fast() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCTGGGCGAGGGGCCAAGGTTGCGCTTTCAGAGCGTTAATATATCCAGTGTAATCAATGATATAATATGGACAGATCATGAAGGTAACTGTTCTTTAG");
        // revcomp
        // CTAAAGAACAGTTACCTTCATGATCTGTCCATATTATATCATTGATTACACTGGATATATTAACGCTCTGAAAGCGCAACCTTGGCCCCTCGCCCAGGGGAGACGAAGGTCCTTCATTACAGTAAATA
        let read = read.as_bytes();
        let mut slice = [0; 128];
        slice.copy_from_slice(read);

        unsafe {
            let register = encode_128_bases_revcomp(&slice);

            let register = reinterpret_as_u8_slice(register);
            assert_eq!(
                register,
                [
                    // divide the string in 4:
                    // CTAAAGAACAGTTACCTTCATGATCTGTCCAT
                    // ATTATATCATTGATTACACTGGATATATTAAC
                    // GCTCTGAAAGCGCAACCTTGGCCCCTCGCCCA
                    // GGGGAGACGAAGGTCCTTCATTACAGTAAATA
                    //
                    // look at the first column: CAGG
                    // then look at the first byte: 0b11110001 <=> GGAC
                    //
                    // each byte represent a column of 4 letters
                    0b11110001, 0b11011010, 0b11101000, 0b11010000, 0b00101000, 0b11110011,
                    0b00001000, 0b01000100, 0b11000001, 0b00111000, 0b00011011, 0b11111110,
                    0b11010010, 0b10001000, 0b01001001, 0b01010001, 0b10010110, 0b10100010,
                    0b01100101, 0b00111000, 0b10111110, 0b10011111, 0b00010000, 0b01011010,
                    0b00010001, 0b11101010, 0b10010011, 0b00111010, 0b00011001, 0b00010001,
                    0b10010000, 0b00000110,
                ]
            );
        }
    }
}

// >>> bin(ord('A'))
// '0b01000001'  # 00 <-> 10 0b10111110 =>
// >>> bin(ord('C'))
// '0b01000011'  # 01 <-> 11 0b10111100 =>
// >>> bin(ord('T'))
// '0b01010100'  # 10 <-> 00 0b10101011 =>
// >>> bin(ord('G'))
// '0b01000111'  # 11 <-> 01 0b10111000 =>

// *

//////
