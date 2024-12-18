use needletail::FastxReader;

pub struct LinesIter {
    data: Box<dyn FastxReader>,
}

impl LinesIter {
    pub fn new(data: Box<dyn FastxReader>) -> Self {
        Self { data }
    }
}

impl Iterator for LinesIter {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.data.next()?.unwrap();
        let sequence = record.seq();
        let sequence_vec = sequence.iter().copied().collect(); // TODO copy
        Some(sequence_vec)
    }
}
