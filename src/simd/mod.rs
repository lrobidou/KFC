// Adapted from Ragnar Groot Koerkamp's minimizers crate to handle ties in a window.
// https://github.com/RagnarGrootKoerkamp/minimizers/tree/master/src/simd

#![allow(dead_code)]
#![cfg_attr(
    not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "avx2"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )),
    deprecated(
        note = "This implementation uses SIMD, make sure you are compiling using `-C target-cpu=native` to get the expected performance."
    )
)]

mod intrinsics;
mod linearize;
mod linearize_backup;
pub mod minimizer;
pub mod nthash;
pub mod packed;
