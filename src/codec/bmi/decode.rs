pub fn decode_32_bases(encoded_bases: u64) -> [u8; 32] {
    let map: [[u8; 2]; 16] = [
        [b'A', b'A'],
        [b'A', b'C'],
        [b'A', b'T'],
        [b'A', b'G'],
        [b'C', b'A'],
        [b'C', b'C'],
        [b'C', b'T'],
        [b'C', b'G'],
        [b'T', b'A'],
        [b'T', b'C'],
        [b'T', b'T'],
        [b'T', b'G'],
        [b'G', b'A'],
        [b'G', b'C'],
        [b'G', b'T'],
        [b'G', b'G'],
    ];
    let mut bytes = [0u8; 32];
    for i in 0..16 {
        let shift = (16 - 1 - i) * 4;
        let two_bases = ((encoded_bases >> shift) & 0xF) as usize;
        let two_bases = unsafe { map.get_unchecked(two_bases) };
        bytes[i * 2..i * 2 + 2].copy_from_slice(two_bases);
    }

    bytes
}

fn decode_up_to_32_bases_slow(encoded_bases: u64, nb_bases: usize) -> Vec<u8> {
    assert!(nb_bases <= 32);
    let mut encoded_bases: u64 = encoded_bases; // >> (2 * (32 - nb_bases));
    let map: [u8; 4] = [b'A', b'C', b'T', b'G'];

    let mut v = vec![];

    for _ in 0..nb_bases {
        let bases_encoded = ((encoded_bases & 0xC000000000000000) >> 62) as usize;
        let ascii_base = map[bases_encoded];
        v.push(ascii_base);
        encoded_bases <<= 2;
    }
    v
}

#[cfg(test)]
pub fn decode_bases(encoded_bases: &[u64], nb_bases: usize) -> Vec<u8> {
    debug_assert_eq!(
        encoded_bases.len(),
        nb_bases / 32 + (((nb_bases % 32) != 0) as usize)
    );
    let nb_iter_full = nb_bases / 32;
    let nb_bases_left_after_full_iter = nb_bases % 32;
    let mut v = vec![];
    encoded_bases
        .iter()
        .take(nb_iter_full)
        .copied()
        .map(decode_32_bases)
        .for_each(|decoded_bases| {
            v.extend_from_slice(&decoded_bases);
        });

    if nb_bases_left_after_full_iter != 0 {
        v.extend_from_slice(&decode_up_to_32_bases_slow(
            encoded_bases[nb_iter_full],
            nb_bases % 32,
        ))
    }
    v
}

pub struct Decoder<'a> {
    nb_bases: usize,
    encoded_bases: &'a [u64],
    // iterates over full blocks of [u8; 32]
    full_values_iterator: std::iter::Take<std::slice::Iter<'a, u64>>,
    full_values_iterator_empty: bool,

    // current full block and position in it
    cache: [u8; 32],
    position_in_cache: usize,

    // rest of bases fater all fulll iterations
    rest: Vec<u8>,
    position_in_rest: usize,
}

impl<'a> Decoder<'a> {
    pub fn new(encoded_bases: &'a [u64], nb_bases: usize) -> Self {
        debug_assert_eq!(
            encoded_bases.len(),
            nb_bases / 32 + (((nb_bases % 32) != 0) as usize)
        );

        let nb_iter_full = nb_bases / 32;
        let mut full_values_iterator: std::iter::Take<std::slice::Iter<u64>> =
            encoded_bases.iter().take(nb_iter_full);

        let first_val = full_values_iterator.next();
        let (cache, full_values_iterator_empty, rest) = if nb_bases == 0 {
            // if there's no base, there is no full u64, and there is no base left after iterating over all the full u64
            let first_full_encoded = [0; 32];
            let full_values_iterator_empty = true;
            let rest_vec = vec![];
            (first_full_encoded, full_values_iterator_empty, rest_vec)
        } else {
            // at least we have some content, but do we have a full u64 ?
            if let Some(x) = first_val {
                // yes => let's decode it
                let first_full_encoded = decode_32_bases(*x);
                let full_values_iterator_empty = false;
                let rest_vec = vec![];
                (first_full_encoded, full_values_iterator_empty, rest_vec)
            } else {
                // no => let's skip it
                let first_full_encoded: [u8; 32] = [0; 32];
                let full_values_iterator_empty = true; // indicates that we skip it
                debug_assert_eq!(encoded_bases.len(), 1); // if there is no full u64, ensure that there there is only one u64 (the partial one)
                let rest_vec = decode_up_to_32_bases_slow(encoded_bases[0], nb_bases);
                (first_full_encoded, full_values_iterator_empty, rest_vec)
            }
        };

        Self {
            nb_bases,
            encoded_bases,
            full_values_iterator,
            full_values_iterator_empty,
            cache,
            position_in_cache: 0,
            rest,
            position_in_rest: 0,
        }
    }
}

impl Iterator for Decoder<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.full_values_iterator_empty {
            let base = self.cache[self.position_in_cache];
            self.position_in_cache += 1;

            // did we reach the end of the cache ?
            if self.position_in_cache >= 32 {
                let next_cache = self.full_values_iterator.next();
                let (cache, full_values_iterator_empty) = if let Some(x) = next_cache {
                    (decode_32_bases(*x), false)
                } else {
                    // note: maybe `self.encoded_bases[self.encoded_bases.len() - 1]` is already decoded
                    // but in this case `self.nb_bases % 32` is 0, so we are not decoding that twice
                    self.rest = decode_up_to_32_bases_slow(
                        self.encoded_bases[self.encoded_bases.len() - 1],
                        self.nb_bases % 32,
                    );
                    ([0; 32], true)
                };
                self.position_in_cache = 0;
                self.cache = cache;
                self.full_values_iterator_empty = full_values_iterator_empty;
            }
            Some(base)
        } else if self.position_in_rest < self.rest.len() {
            let base = self.rest[self.position_in_rest];
            self.position_in_rest += 1;
            Some(base)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_decode_32_bases() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCC");
        let read = read.as_bytes();
        let encoded = 0b10001010_10000110_11100000_10110000_11110001_01101001_11100110_01010101;
        let decoded = decode_32_bases(encoded);
        assert_eq!(read, decoded);
    }

    #[test]
    fn test_decode_33_bases() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCT");
        let read = read.as_bytes();
        let encoded = vec![
            0b10001010_10000110_11100000_10110000_11110001_01101001_11100110_01010101,
            0b10000000_00000000_00000000_00000000_00000000_00000000_00000000_00000000,
        ];
        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read, decoded);
    }

    #[test]
    fn test_decode_29_bases() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTC");
        let read = read.as_bytes();
        let encoded =
            vec![0b10001010_10000110_11100000_10110000_11110001_01101001_11100110_01000000];
        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read, decoded);
    }

    #[test]
    fn test_decode_empty() {
        let read = String::from("");
        let read = read.as_bytes();
        let encoded = vec![];
        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read, decoded);
    }
}
