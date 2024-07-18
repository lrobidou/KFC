# KFC

## TODO:
- [x] using &str instead of String to prevent copies
- [x] using Vec<u8> instead of `String` in the data structures
- [x] implement search
- [] implement streaming search
- [] use kff as possible output
- [] implement kmer iterator
- [x] streaming of superkmer

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