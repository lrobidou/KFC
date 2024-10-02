use serde::de::{SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::marker::PhantomData;
use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::{cmp, fmt};

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
            let chunks = self.chunks();
            for chunk in <[Arc<RwLock<T>>; NB_BUCKETS] as Clone>::clone(chunks).into_iter() {
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
        let data: [Arc<RwLock<T>>; NB_BUCKETS] =
            std::array::from_fn(|_i| Arc::new(RwLock::new(function())));
        Self { data }
    }

    pub fn get_from_id_usize(&self, id: usize) -> Arc<RwLock<T>> {
        let id = id % NB_BUCKETS;
        self.data[id].clone()
    }

    pub fn get_from_id_u64(&self, id: u64) -> Arc<RwLock<T>> {
        let id = usize::try_from(id).unwrap();
        self.get_from_id_usize(id)
    }

    pub fn chunks(&self) -> &[Arc<RwLock<T>>; NB_BUCKETS] {
        &self.data
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn acquire_write_locks(
        &self,
        idx0: u64, // current minimizer will be locked for sure
        idx1: u64, // previous minimizer
        idx2: u64, // next minimizer
    ) -> (
        RwLockWriteGuard<T>,
        Option<RwLockReadGuard<T>>,
        Option<RwLockReadGuard<T>>,
        LockPosition,
        LockPosition,
    ) {
        let idx0 = usize::try_from(idx0).unwrap() % NB_BUCKETS;
        let idx1 = usize::try_from(idx1).unwrap() % NB_BUCKETS;
        let idx2 = usize::try_from(idx2).unwrap() % NB_BUCKETS;

        let mut lock_0 = None;
        let mut lock_1 = None;
        let mut lock_2 = None;

        let mut lock_1_position = LockPosition::ThisLock;
        let mut lock_2_position = LockPosition::ThisLock;

        // Create a list of indices with their original position
        let mut indices = [(idx0, 0), (idx1, 1), (idx2, 2)];

        // Sort by index to avoid deadlocks
        indices.sort_by_key(|&(idx, _)| idx);

        // Iterate over the sorted indices and acquire the necessary locks
        // let mut previous_val_and_pos = None;
        for (idx, original_pos) in indices.into_iter() {
            if idx == idx0 {
                if lock_0.is_none() {
                    lock_0 = Some(self.data[idx].write().unwrap());
                }
                if original_pos == 0 {
                    // nothing to do
                } else if original_pos == 1 {
                    lock_1_position = LockPosition::CurrentLock;
                } else {
                    lock_2_position = LockPosition::CurrentLock;
                }
            } else if idx == idx1 {
                if lock_1.is_none() {
                    // lock_id_1 = Some(original_pos);
                    lock_1 = Some(self.data[idx].read().unwrap());
                }
                if original_pos == 0 {
                    // TODO please approve me
                    // impossible
                    panic!()
                } else if original_pos == 1 {
                    lock_1_position = LockPosition::ThisLock;
                } else {
                    lock_2_position = LockPosition::OtherLock;
                }
            } else {
                assert!(idx == idx2);
                assert!(lock_2.is_none());
                if lock_2.is_none() {
                    // lock_id_2 = Some(original_pos);
                    lock_2 = Some(self.data[idx].read().unwrap());
                }
                if original_pos == 0 {
                    panic!()
                    // TODO please approve me
                    // impossible
                } else if original_pos == 1 {
                    panic!()
                    // TODO please approve me
                    // impossible
                } else {
                    lock_2_position = LockPosition::ThisLock;
                }
            }
        }

        (
            lock_0.unwrap(),
            lock_1,
            lock_2,
            lock_1_position,
            lock_2_position,
        )
    }

    /// Acquire two locks in read mode. The first one is always acquired.
    /// If the second one is equal to the same one, it is set to None.
    /// The locks are sorted before TODO acquiring them.
    /// This prevents deadlocks.
    pub fn acquire_two_locks_read_mode(
        &self,
        idx0: usize, // will be locked for sure
        idx1: usize, // might not be locked
    ) -> (RwLockReadGuard<T>, Option<RwLockReadGuard<T>>) {
        let idx0 = idx0 % NB_BUCKETS;
        let idx1 = idx1 % NB_BUCKETS;

        match idx0.cmp(&idx1) {
            cmp::Ordering::Less => {
                let lock_1 = self.data[idx0].read().unwrap();
                let lock_2 = Some(self.data[idx1].read().unwrap());
                (lock_1, lock_2)
            }
            cmp::Ordering::Equal => {
                let lock_1 = self.data[idx0].read().unwrap();
                let lock_2 = None;
                (lock_1, lock_2)
            }
            cmp::Ordering::Greater => {
                let lock_2 = Some(self.data[idx1].read().unwrap());
                let lock_1 = self.data[idx0].read().unwrap();
                (lock_1, lock_2)
            }
        }
    }
}

#[derive(Debug)]
pub enum LockPosition {
    ThisLock,
    CurrentLock,
    OtherLock,
}
