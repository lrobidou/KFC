mod decode;
mod encode;
mod revcomp;

pub use decode::Decoder;
pub use encode::{Encoder, RevCompEncoder};
pub use revcomp::{FusedReverseIterator, RevCompIter};

use revcomp::revcomp_32_bases;

#[cfg(test)]
mod tests {
    use decode::{decode_32_bases, decode_bases};
    use encode::{encode_32_bases_revcomp, encode_bases};
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_encode_32_bases_revcomp() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCC");
        let read_revcomp = String::from("GGGGAGACGAAGGTCCTTCATTACAGTAAATA");
        let read = read.as_bytes();
        let read_revcomp = read_revcomp.as_bytes();

        let encoded = encode_32_bases_revcomp(read.try_into().unwrap());
        let decoded = decode_32_bases(encoded);
        assert_eq!(decoded, read_revcomp);
    }

    #[test]
    fn test_encode_bases_32_bases() {
        let read = String::from("TTGGGAGAGTGTCTGAACGTGCGCTAGAGTGGTGACTCTAAGTGCAATAACCCTTCCTTTGGGT");
        let read = read.as_bytes();

        let encoded = encode_bases(read);
        let decoded = decode_bases(&encoded, read.len());
        assert_eq!(decoded, read);
    }

    #[test]
    fn test_encode_bases() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTAAGGACCTTCGTCTCCCCGG");
        let read = read.as_bytes();

        let encoded = encode_bases(read);
        let decoded = decode_bases(&encoded, read.len());
        assert_eq!(decoded, read);
    }

    #[test]
    fn test_decoder() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTAAGGACCTTCGTCTCCCCGG");
        let read = read.as_bytes();

        let encoded = encode_bases(read);

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read, decoded);
    }

    #[test]
    fn test_encode_decode() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTAAGGACCTTCGTCTCCCCGG");
        let read = read.as_bytes();

        let encoder = Encoder::new(read);
        let encoded = encoder.collect_vec();

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read, decoded);
    }

    #[test]
    fn test_revcomp_encoder() {
        let read = String::from("TATTTACTGTAATGAAGGACCTTCGTCTCCCCGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTGGGGAGACGAAGGTCCTTAAGGACCTTCGTCTCCCCGG");
        let read = read.as_bytes();
        let read_revcomp = String::from("CCGGGGAGACGAAGGTCCTTAAGGACCTTCGTCTCCCCAAGGACCTTCGTCTCCCCAAGGACCTTCGTCTCCCCGGGGAGACGAAGGTCCTTCATTACAGTAAATA");
        let read_revcomp = read_revcomp.as_bytes();

        let encoder = RevCompEncoder::new(read);
        let encoded = encoder.collect_vec();

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read_revcomp, decoded);
    }

    #[test]
    fn test_revcomp_encoder_32_bases() {
        let read = String::from("GCGCAACCTATTTAGGTTTTCCTAGAGGTCGA");
        let read = read.as_bytes();
        let read_revcomp = String::from("TCGACCTCTAGGAAAACCTAAATAGGTTGCGC");
        let read_revcomp = read_revcomp.as_bytes();

        let encoder = RevCompEncoder::new(read);
        let encoded = encoder.collect_vec();

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read_revcomp, decoded);
    }

    #[test]
    fn test_revcomp_encoder_10_bases() {
        let read = String::from("GTCTGGACTA");
        let read = read.as_bytes();
        let read_revcomp = String::from("TAGTCCAGAC");
        let read_revcomp = read_revcomp.as_bytes();

        let encoder = RevCompEncoder::new(read);
        let encoded = encoder.collect_vec();

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read_revcomp, decoded);
    }

    #[test]
    fn test_revcomp_encoder_0_bases() {
        let read = String::from("");
        let read = read.as_bytes();
        let read_revcomp = String::from("");
        let read_revcomp = read_revcomp.as_bytes();

        let encoder = RevCompEncoder::new(read);
        let encoded = encoder.collect_vec();

        let decoder = Decoder::new(&encoded, read.len());
        let decoded = decoder.collect_vec();
        assert_eq!(read_revcomp, decoded);
    }

    #[test]
    fn revcomp_encoder() {
        let read = String::from("GTTCCAT");
        let revcomp_read = String::from("ATGGAAC");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            assert_eq!(
                encoded,
                //  G T T C C A T
                [0b1110100101001000000000000000000000000000000000000000000000000000]
            );

            let revcomp = RevCompEncoder::new(bytes).collect_vec();
            assert_eq!(
                revcomp,
                //  T A G G A A C
                [0b0010111100000100000000000000000000000000000000000000000000000000]
            );
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }

    #[test]
    fn revcomp_encoder_2() {
        let read = String::from("ATCGGCGCATCG");
        let revcomp_read = String::from("CGATGCGCCGAT");
        let bytes = read.as_bytes();
        {
            let encoder = Encoder::new(bytes);
            let encoded = encoder.collect_vec();
            assert_eq!(
                encoded,
                //  A T C G G C G C A T C G
                [0b0010011111011101001001110000000000000000000000000000000000000000]
            );

            let revcomp = RevCompEncoder::new(bytes).collect_vec();
            assert_eq!(
                revcomp,
                //  C G A T G C G C C G A T
                [0b0111001011011101011100100000000000000000000000000000000000000000]
            );
            let decoded = Decoder::new(&revcomp, bytes.len()).collect_vec();
            assert_eq!(revcomp_read.as_bytes(), decoded);
        }
    }
}
