use crate::Minimizer;
use serde::de::{SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::marker::PhantomData;
use std::sync::{Arc, RwLock};

pub const NB_PARALLEL_CHUNK: usize = 255;

pub struct Parallel<T: PartialEq + Serialize + for<'a> Deserialize<'a>> {
    data: [Arc<RwLock<T>>; NB_PARALLEL_CHUNK],
}

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Clone for Parallel<T> {
    fn clone(&self) -> Parallel<T> {
        Parallel {
            data: self.data.clone(),
        }
    }
}

// Implement Serialize for Parallel
impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Serialize for Parallel<T> {
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

// Implement Deserialize for Parallel
impl<'de, T> Deserialize<'de> for Parallel<T>
where
    T: PartialEq + Serialize + for<'a> Deserialize<'a>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct ParallelVisitor<T>(PhantomData<T>);

        impl<'de, T> Visitor<'de> for ParallelVisitor<T>
        where
            T: PartialEq + Serialize + for<'a> Deserialize<'a>,
        {
            type Value = Parallel<T>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("an array of length NB_PARALLEL_CHUNK")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let data: [Arc<RwLock<T>>; NB_PARALLEL_CHUNK] = std::array::from_fn(|_| {
                    Arc::new(RwLock::new(seq.next_element().unwrap().unwrap()))
                });

                Ok(Parallel { data })
            }
        }

        deserializer.deserialize_seq(ParallelVisitor(PhantomData))
    }
}

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> PartialEq for Parallel<T> {
    fn eq(&self, other: &Parallel<T>) -> bool {
        for i in 0..NB_PARALLEL_CHUNK {
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

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Parallel<T> {
    pub fn new<F>(function: F) -> Self
    where
        F: Fn() -> T,
    {
        let data = std::array::from_fn(|_i| Arc::new(RwLock::new(function())));
        Self { data }
    }

    pub fn get_from_minimizer(&self, minimizer: Minimizer) -> Arc<RwLock<T>> {
        let index = minimizer % Minimizer::try_from(NB_PARALLEL_CHUNK).unwrap();
        let index = usize::try_from(index).unwrap();
        // TODO unchecked
        self.data[index].clone()
    }

    pub fn chunks(&self) -> &[Arc<RwLock<T>>; NB_PARALLEL_CHUNK] {
        &self.data
    }
}
