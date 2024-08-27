## TODO:
- [x] using &str instead of String to prevent copies
- [x] using Vec<u8> instead of `String` in the data structures
- [x] implement search
- [ ] implement streaming search
- [x] use kff as possible output
- [x] implement k-mer iterator
- [x] streaming of k-mers
- [ ] optimize iteration of minimizer for streaming of k-mers
- [ ] discuss what to do when a single k-mer is indexed: there is no previous nor next k-mer
- [ ] debug: we miss the first and last superkmer if t = 1 
- [ ] use version of kff from the crate instead of my own 
- [ ] make it possible to do a single pass 
- [ ] allow non canonical k-mers (?)
- [ ] reference the paper in the README
- [ ] test large value of k
- [ ] make a paper branch (?)
- [ ] test nightly
- [ ] self mut on MashMap::get_mut_iter ?
- [ ] better parallelization than locking the whole minimizer bucket ?
- [ ] do we clone Buckets or do we use Arc ?

## bug lefts
We miss the first and last superkmers:
```bash
cargo run -- build -k 99 -m 21 -t 1 --input data/U00096.3.fasta --output index_64.kfc --check-kmc data/99mers.txt
cargo run -- dump --input-index index_64.kfc --output-text index_64.txt
# sort the two files
sort index_64.txt > dump_sorted.txt 
sort 99mers.txt > 99_sort.txt
diff dump_sorted.txt 99_sort.txt  # 74 k-mers missing from the output of KFC
```
Alternatively:
```bash
cargo run -- build -k 99 -m 21 -t 1 --input data/U00096.3.fasta --output index_64.kfc --check-kmc data/99mers.txt
cargo run -- dump --input-index index_64.kfc --output-kff index.kff
cargo run -- kff-dump --input-kff index.kff --output-text kff_text.txt

# sort the two files
sort kff_text.txt > kff_sorted.txt
sort 99mers.txt > 99_sort.txt
diff kff_sorted.txt 99_sort.txt  # 74 k-mers missing from the output of KFC
```
# KFC

![tests](https://github.com/lrobidou/KFC/workflows/tests/badge.svg)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/lrobidou/KFC//blob/main/LICENSE)

KFC is an implementation of a canonical k-mer counter using a k-mer representation introduced in XXX the paper XXX.

This representation allows to:
- use a (very, very) large k value (as big as 2^(regirster size)) XXX more tests on big values XXX
- be space efficient (output index is typically ~x XXX x XXX times smaller than KMC)
- be fast (indexation is typically ~x XXX x XXX faster than KMC)
KFC can filter k-mer based on their count and only index the k-mers above that count.

## If you are a reviewer
We created a branch `paper` that will stay consistant with the paper.

Please checkout to the `paper` branch:
```bash
git clone https://github.com/lrobidou/KFC
git checkout paper XXX do the branch XXX
cd KFC
cargo build --release
```

## Install

First, [install rust](https://www.rust-lang.org/learn/get-started).

Then install KFC:
```bash
git clone https://github.com/lrobidou/KFC
cd KFC
cargo build --release
```
XXX install in path ? XXX

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