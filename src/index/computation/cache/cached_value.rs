#[derive(Debug, Clone)]
pub struct CachedValue {
    id_bucket: usize,
    id_hk: usize,
    is_large: bool,
}

impl CachedValue {
    pub fn new(id_bucket: usize, id_hk: usize, is_large: bool) -> Self {
        Self {
            id_bucket,
            id_hk,
            is_large,
        }
    }

    pub fn get_id_bucket(&self) -> usize {
        self.id_bucket
    }

    pub fn get_id_hk(&self) -> usize {
        self.id_hk
    }

    pub fn get_is_large(&self) -> bool {
        self.is_large
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_cache() {
        let id_bucket: usize = 5;
        let id_hk: usize = 9;
        let is_large: bool = false;
        let cache_value = CachedValue::new(id_bucket, id_hk, is_large);

        assert_eq!(cache_value.get_id_bucket(), id_bucket);
        assert_eq!(cache_value.get_id_hk(), id_hk);
        assert_eq!(cache_value.get_is_large(), is_large);
    }
}
