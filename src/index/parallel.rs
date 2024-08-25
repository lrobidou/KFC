use crate::Minimizer;
use serde::de::{SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::marker::PhantomData;
use std::sync::{Arc, RwLock};

pub const NB_PARALELL_CHUNK: usize = 255;

pub struct Paralell<T: PartialEq + Serialize + for<'a> Deserialize<'a>> {
    data: [Arc<RwLock<T>>; NB_PARALELL_CHUNK],
}

// Implement Serialize for Paralell
impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Serialize for Paralell<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        {
            let mut seq = serializer.serialize_seq(Some(self.chunks().len()))?;
            for chunk in self.chunks() {
                let chunk = chunk.read().unwrap();
                seq.serialize_element(&*chunk)?;
            }
            seq.end()
        }
    }
}

// Implement Deserialize for Paralell
impl<'de, T> Deserialize<'de> for Paralell<T>
where
    T: PartialEq + Serialize + for<'a> Deserialize<'a>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ParalellVisitor<T>(PhantomData<T>);

        impl<'de, T> Visitor<'de> for ParalellVisitor<T>
        where
            T: PartialEq + Serialize + for<'a> Deserialize<'a>,
        {
            type Value = Paralell<T>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("an array of length NB_PARALELL_CHUNK")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let data: [Arc<RwLock<T>>; NB_PARALELL_CHUNK] = std::array::from_fn(|_| {
                    Arc::new(RwLock::new(seq.next_element().unwrap().unwrap()))
                });

                Ok(Paralell { data })
            }
        }

        deserializer.deserialize_seq(ParalellVisitor(PhantomData))
    }
}

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> PartialEq for Paralell<T> {
    fn eq(&self, other: &Paralell<T>) -> bool {
        for i in 0..NB_PARALELL_CHUNK {
            // TODO get_unchecked
            let data = self.data[i].read().unwrap();
            let other_data = other.data[i].read().unwrap();
            if *data != *other_data {
                return false;
            }
        }
        true
    }
}

// pub mod mac {
//     macro_rules! read_lock {
//         ($var:expr, $minimizer:expr) => {{
//             let chunk = $var.get_from_minimizer($minimizer);
//             let lock = chunk.read().unwrap();
//             lock
//         }};
//     }

//     macro_rules! write_lock {
//         ($var:expr, $minimizer:expr) => {{
//             $var.get_from_minimizer($minimizer).write().unwrap()
//         }};
//     }

//     pub(crate) use read_lock;
//     pub(crate) use write_lock;
// }

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Paralell<T> {
    pub fn new(function: fn() -> T) -> Self {
        let data = std::array::from_fn(|_i| Arc::new(RwLock::new(function())));
        Self { data }
    }

    pub fn get_from_minimizer(&self, minimizer: Minimizer) -> Arc<RwLock<T>> {
        let index = minimizer % Minimizer::try_from(NB_PARALELL_CHUNK).unwrap();
        let index = usize::try_from(index).unwrap();
        // TODO unchecked
        self.data[index].clone()
    }

    pub fn chunks(&self) -> &[Arc<RwLock<T>>; NB_PARALELL_CHUNK] {
        &self.data
    }
}
