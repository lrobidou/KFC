mod extended_hyperkmers;
mod hyperkmers_counts;
mod superkmers_count;

pub use extended_hyperkmers::ExtendedHyperkmers;
pub use hyperkmers_counts::{search_exact_hyperkmer_match, HKCount, HKMetadata};
pub use superkmers_count::SuperKmerCounts;
