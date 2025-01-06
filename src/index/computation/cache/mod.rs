mod cache_misere;
mod cached_value;

pub use cache_misere::CacheMisere;
pub use cached_value::CachedValue;

// #[macro_export]
// macro_rules! to_cache {
//     ($v:ident) => {{
//         $v.clear();
//         $v.into_iter().map(|_| unreachable!()).collect_vec().into()
//     }};
// }
