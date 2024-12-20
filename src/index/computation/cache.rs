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
