use crate::Minimizer;
use serde::de::{SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::marker::PhantomData;
use std::sync::{Arc, RwLock};

pub const NB_BUCKETS: usize = 255;

pub struct Buckets<T: PartialEq + Serialize + for<'a> Deserialize<'a>> {
    data: [Arc<RwLock<T>>; NB_BUCKETS],
}

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Clone for Buckets<T> {
    fn clone(&self) -> Buckets<T> {
        Buckets {
            data: self.data.clone(),
        }
    }
}

// Implement Serialize for Buckets
impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Serialize for Buckets<T> {
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

// Implement Deserialize for Buckets
impl<'de, T> Deserialize<'de> for Buckets<T>
where
    T: PartialEq + Serialize + for<'a> Deserialize<'a>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct BucketsVisitor<T>(PhantomData<T>);

        impl<'de, T> Visitor<'de> for BucketsVisitor<T>
        where
            T: PartialEq + Serialize + for<'a> Deserialize<'a>,
        {
            type Value = Buckets<T>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("an array of length NB_BUCKETS")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let data: [Arc<RwLock<T>>; NB_BUCKETS] = std::array::from_fn(|_| {
                    Arc::new(RwLock::new(seq.next_element().unwrap().unwrap()))
                });

                Ok(Buckets { data })
            }
        }

        deserializer.deserialize_seq(BucketsVisitor(PhantomData))
    }
}

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> PartialEq for Buckets<T> {
    fn eq(&self, other: &Buckets<T>) -> bool {
        for i in 0..NB_BUCKETS {
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

impl<T: PartialEq + Serialize + for<'a> Deserialize<'a>> Buckets<T> {
    pub fn new<F>(function: F) -> Self
    where
        F: Fn() -> T,
    {
        let data = std::array::from_fn(|_i| Arc::new(RwLock::new(function())));
        Self { data }
    }

    pub fn get_from_minimizer(&self, minimizer: Minimizer) -> Arc<RwLock<T>> {
        let index = minimizer % Minimizer::try_from(NB_BUCKETS).unwrap();
        let index = usize::try_from(index).unwrap();
        // TODO unchecked
        self.data[index].clone()
    }

    pub fn chunks(&self) -> &[Arc<RwLock<T>>; NB_BUCKETS] {
        &self.data
    }
}
