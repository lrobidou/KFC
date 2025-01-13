use crate::subsequence::{NoBitPacked, Subsequence};

pub type AtypicalHyperkmerPart = (usize, Vec<u64>);

use serde::{
    de::{SeqAccess, Visitor},
    Deserialize, Deserializer, Serialize, Serializer,
};

#[derive(PartialEq)]
pub struct AtypicalHyperkmerParts {
    data: boxcar::Vec<AtypicalHyperkmerPart>,
}

impl Serialize for AtypicalHyperkmerParts {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeSeq;

        let mut seq = serializer.serialize_seq(Some(self.data.count()))?;
        for (_index, array) in &self.data {
            seq.serialize_element(array)?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for AtypicalHyperkmerParts {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "snake_case")]
        enum Field {
            Data,
        }

        struct AtypicalHyperkmerPartsVisitor;

        impl<'de> Visitor<'de> for AtypicalHyperkmerPartsVisitor {
            type Value = AtypicalHyperkmerParts;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                formatter.write_str("hyperkmer parts bucket")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<AtypicalHyperkmerParts, V::Error>
            where
                V: SeqAccess<'de>,
            {
                let arrays: Vec<(usize, Vec<u64>)> = seq
                    .next_element()?
                    .ok_or_else(|| serde::de::Error::invalid_length(4, &self))?;
                let data = arrays.into_iter().collect();
                Ok(AtypicalHyperkmerParts { data })
            }
        }

        const FIELDS: &[&str] = &["data"];
        deserializer.deserialize_struct(
            "HyperkmerPartsBucket", // TODO name ?
            FIELDS,
            AtypicalHyperkmerPartsVisitor,
        )
    }
}

impl AtypicalHyperkmerParts {
    pub fn new() -> Self {
        Self {
            data: boxcar::Vec::new(),
        }
    }

    pub fn get(&self, index: usize) -> Option<&(usize, Vec<u64>)> {
        self.data.get(index)
    }

    // TODO make an iterator instead
    pub fn get_data(&self) -> &boxcar::Vec<(usize, Vec<u64>)> {
        &self.data
    }

    /// Adds `new_ext_hyperkmer` into `large_hyperkmers`.
    ///
    /// Returns the position of `new_ext_hyperkmer` in `large_hyperkmers`
    pub fn add_new_large_hyperkmer(&self, new_ext_hyperkmer: &Subsequence<NoBitPacked>) -> usize {
        let nb_base = new_ext_hyperkmer.len();

        let size_needed_for_slice = nb_base / 32 + (nb_base % 32 != 0) as usize;

        let mut dest_slice = vec![0; size_needed_for_slice];
        new_ext_hyperkmer
            .to_canonical()
            .dump_as_2bits(&mut dest_slice);

        self.data.push((nb_base, dest_slice))
    }
}

// #[cfg(any(debug_assertions, test))]
impl AtypicalHyperkmerParts {
    pub fn len(&self) -> usize {
        self.data.count()
    }
}
