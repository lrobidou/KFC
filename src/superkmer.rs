use super::superkmers_computation::is_canonical;
use xxhash_rust::const_xxh3::xxh3_64;
pub fn reverse_complement<'a>(seq: &'a str) -> String {
    seq.as_bytes()
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            _ => *base as char,
        })
        .collect()
}

// TODO is there no copy here ? What is the cost of moving references ?
pub fn reverse_complement_no_copy<'a>(
    seq: &'a str,
) -> std::iter::Map<std::iter::Rev<std::slice::Iter<'_, u8>>, fn(&u8) -> u8> {
    seq.as_bytes().iter().rev().map(|base| match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => *base,
    })
}

// TODO duplication of `read` and `same_orientation` when used in Superkmer
/// Represents a subsequence, possibly in reverse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SubsequenceMetadata<'a> {
    read: &'a str,
    start: usize,
    end: usize,
    same_orientation: bool,
}

impl<'a> SubsequenceMetadata<'a> {
    pub fn new(read: &'a str, start: usize, end: usize, same_orientation: bool) -> Self {
        assert!(start <= read.len());
        assert!(end <= read.len());
        Self {
            read,
            start,
            end,
            same_orientation,
        }
    }
    pub fn whole_string(read: &'a str) -> Self {
        Self::new(read, 0, read.len(), true)
    }

    pub fn change_orientation(&self) -> Self {
        Self {
            read: self.read,
            start: self.start,
            end: self.end,
            same_orientation: !self.same_orientation,
        }
    }

    pub fn change_orientation_if(&self, cond: bool) -> Self {
        if cond {
            Self {
                read: self.read,
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
            }
        } else {
            Self {
                read: self.read,
                start: self.start,
                end: self.end,
                same_orientation: self.same_orientation,
            }
        }
    }

    pub fn start(&self) -> usize {
        self.start
    }
    pub fn end(&self) -> usize {
        self.end
    }
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    // TODO there is a copy here
    pub fn is_canonical(&self) -> bool {
        let subsequence = self.to_string();
        is_canonical(&subsequence)
    }

    pub fn to_canonical(&self) -> Self {
        let subsequence = &self.read[self.start..self.end];

        if is_canonical(subsequence) == self.same_orientation {
            *self
        } else {
            Self {
                read: self.read,
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
            }
        }
    }

    // TODO remove
    pub fn to_string(&self) -> String {
        let subsequence = &self.read[self.start..self.end];
        if self.same_orientation {
            subsequence.into()
        } else {
            reverse_complement(subsequence)
        }
    }

    // TODO remove
    pub fn to_canonical_string(&self) -> String {
        let subsequence = &self.read[self.start..self.end];
        if is_canonical(subsequence) {
            subsequence.into()
        } else {
            reverse_complement(subsequence)
        }
    }

    // TODO remove copy made by `reverse_complement`
    pub fn equal_str(&self, other: &str) -> bool {
        self.to_string() == other
    }

    // TODO remove copy here
    pub fn ends_with(&self, other: &SubsequenceMetadata) -> bool {
        let x = self.to_string();
        let y = other.to_string();
        x.ends_with(&y)
    }

    // TODO remove copy here
    pub fn starts_with(&self, other: &SubsequenceMetadata) -> bool {
        let x = self.to_string();
        let y = other.to_string();
        x.starts_with(&y)
    }

    /// Extract a subsequence
    /// Equivalent to:
    /// ```
    /// Self::whole_string(self.to_string()[start..end])
    /// ```
    pub fn subsequence(&self, start: usize, end: usize) -> SubsequenceMetadata {
        // println!("taking a subsequence of {}", self.to_string());
        // println!("{:?}", (start, end));
        if self.same_orientation {
            Self {
                read: self.read,
                start: self.start + start,
                end: self.start + end,
                same_orientation: self.same_orientation,
            }
        } else {
            Self {
                read: self.read,
                start: self.end - end,
                end: self.end - start,
                same_orientation: self.same_orientation,
            }
        }
    }

    pub fn common_prefix_length(&self, other: &SubsequenceMetadata) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.read[self.start..self.end].bytes();
            let mut y = other.read[other.start..other.end].bytes();
            iter_prefix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.read[self.start..self.end].bytes();
            let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
            iter_prefix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
            let mut y = other.read[other.start..other.end].bytes();
            iter_prefix_len(&mut x, &mut y)
        } else {
            //f !self.same_orientation && !other.same_orientation {
            let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
            let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
            iter_prefix_len(&mut x, &mut y)
        }
    }

    pub fn common_suffix_length(&self, other: &SubsequenceMetadata) -> usize {
        if self.same_orientation && other.same_orientation {
            let mut x = self.read[self.start..self.end].bytes();
            let mut y = other.read[other.start..other.end].bytes();
            iter_suffix_len(&mut x, &mut y)
        } else if self.same_orientation && !other.same_orientation {
            let mut x = self.read[self.start..self.end].bytes();
            let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
            iter_suffix_len(&mut x, &mut y)
        } else if !self.same_orientation && other.same_orientation {
            let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
            let mut y = other.read[other.start..other.end].bytes();
            iter_suffix_len(&mut x, &mut y)
        } else {
            //f !self.same_orientation && !other.same_orientation {
            let mut x = reverse_complement_no_copy(&self.read[self.start..self.end]);
            let mut y = reverse_complement_no_copy(&other.read[other.start..other.end]);
            iter_suffix_len(&mut x, &mut y)
        }
    }
}

fn iter_prefix_len(x: &mut impl Iterator<Item = u8>, y: &mut impl Iterator<Item = u8>) -> usize {
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}

fn iter_suffix_len(
    x: &mut impl DoubleEndedIterator<Item = u8>,
    y: &mut impl DoubleEndedIterator<Item = u8>,
) -> usize {
    let mut x = x.rev().into_iter();
    let mut y = y.rev().into_iter();
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}
#[derive(Debug, Clone, Copy, PartialEq)]

pub struct Superkmer<'a> {
    pub read: &'a str,
    pub minimizer: SubsequenceMetadata<'a>,
    pub superkmer: SubsequenceMetadata<'a>,
}

impl<'a> Superkmer<'a> {
    pub fn new(
        read: &'a str,
        start_mini: usize,
        end_mini: usize,
        start_sk: usize,
        end_sk: usize,
        same_orientation: bool,
    ) -> Self {
        Self {
            read,
            minimizer: SubsequenceMetadata {
                read,
                start: start_mini,
                end: end_mini,
                same_orientation,
            },
            superkmer: SubsequenceMetadata {
                read,
                start: start_sk,
                end: end_sk,
                same_orientation,
            },
        }
    }
    pub fn hash_superkmer(&self) -> u64 {
        xxh3_64(self.superkmer.to_string().as_bytes())
    }

    // TODO not String anymore
    pub fn get_minimizer(&self) -> String {
        self.minimizer.to_string()
    }

    pub fn m(&self) -> usize {
        self.minimizer.len()
    }

    pub fn is_canonical_in_the_read(&self) -> bool {
        self.minimizer.same_orientation
    }

    pub fn print(&self) -> String {
        format!(
            "{} ({} {})",
            &self.superkmer.read[self.superkmer.start..self.superkmer.end],
            self.minimizer.to_string(),
            self.minimizer.same_orientation
        )
    }
}
