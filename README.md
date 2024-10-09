# KFC

![tests](https://github.com/lrobidou/KFC/workflows/tests/badge.svg)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/lrobidou/KFC//blob/main/LICENSE)

KFC is an implementation of a canonical k-mer counter using a k-mer representation introduced in XXX the paper XXX.

This representation allows to:
- use a (very, very) large k value (as big as 2^(regirster size)) XXX more tests on big values XXX
- be space efficient (output index is typically ~x XXX x XXX times smaller than KMC)
- be fast (indexation is typically ~x XXX x XXX faster than KMC)
KFC can filter k-mer based on their count and only index the k-mers above that count.

## Install

If you have not installed Rust yet, please visit [rustup.rs](https://rustup.rs/) to install it.

Then clone this repository and build KFC using:
```sh
git clone https://github.com/lrobidou/KFC
cd KFC
RUSTFLAGS="-C target-cpu=native" cargo +nightly build --release -F nightly
```
Make sure to set `RUSTFLAGS="-C target-cpu=native"` to use the fastest instructions available on your architecture.

If you cannot use Rust nightly, you can also build KFC in stable mode (which may be slightly slower):
```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release
```
This will create a binary located at `target/release/kfc`.

## Run
### Build a KFC index
The first step to any KFC usage is to build a KFC index.
```bash
cargo run --release -- build -k <k>> -m <m> -t <threshold_count> --input <file>.fasta --output <index>.kfc
```

### Dump a KFC index to text
Once the KFC index is computed, it is possible to dump it to text. The k-mers are *not* ordered.
```bash
cargo run --release -- dump --input-index <index>.kfc --output-text <kmers.txt>
```

### Dump a KFC index to the k-mer file format (KFF)
KFC supports the k-mer file format (see [Dufresne et al, The K-mer File Format: a standardized and compact disk representation of sets of k-mers](https://doi.org/10.1093/bioinformatics/btac528)).
As such, it is possible to dump a KFC index into a KFF file.
The count of each k-mer is encoded in the KFF file.
```bash
cargo run --release -- dump --input-index <index>.kfc --output-kff <index>.kff
```

### Dump a KFF to text
**Warning:** KFC only handles KFF files built by KFC.

Reading the KFF file produced by KFC should be possible with any implementation supporting KFF, but we recommand relying on KFC for this task. Indeed, a KFF built by KFC respects some assumptions on the count of k-mers, which can be used to dump the KFF file with a lower memory consumption. This also means that files not respecting these assumptions would produce invalid count if dumped by KFC.

```bash
cargo run --release -- kff-dump --input-kff <index>.kff --output-text <index>.txt
```

# Tests
# data sources
Place [https://www.ebi.ac.uk/ena/browser/view/U00096](https://www.ebi.ac.uk/ena/browser/view/U00096) in data/U00096.3.fasta

# Contributing

## Testing

### Install dependencies

* [cargo-nextest](https://nexte.st/)
* [cargo-llvm-cov](https://crates.io/crates/cargo-llvm-cov)
```bash
cargo install cargo-llvm-cov
cargo install cargo-nextest
```

### Generate report on terminal

At the workspace root:
```bash
cargo nextest run
```

### Generating coverage

See at the beginning of the bellow script the packages to install.

At the workspace root:

```bash
./coverage.sh
cargo llvm-cov --open # open HTML report into the navigator
```

It also generates a `.lcov.info` lcov file.

## If you have a performance issue:
### Flamegraph
You can generate a flamegraph using [cargo flamegraph](https://github.com/flamegraph-rs/flamegraph):
```bash
sudo apt install -y linux-perf  # for debian distributions
cargo install flamegraph
echo -1 | sudo tee /proc/sys/kernel/perf_event_paranoid  # do this each time you reboot
cargo flamegraph --image-width 30000 -- build -k 99 -m 21 -t 1 --input data/U00096.3.fasta # change the width to suit your need
```
This graph allows you to quickly check wich function is taking the more time: the longer a function takes in total, the longest the associated rectangle.

### Using vtune
If you you accept to run arbitrary closed-source binary on you system, you can install vtune from here : https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler-download.html?operatingsystem=linux&linux-install-type=offline. You do not need to create an account, just continue as guest.
Execute the script, it will install a binary. Locate it and execture it.
```bash
~/intel/oneapi/vtune/2024.2/bin64/vtune-gui
```
vtune is able to give much more informaton, including a nicer version of the flamegraph.