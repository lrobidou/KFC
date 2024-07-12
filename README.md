# KFC

## TODO:
- [x] using &str instead of String to prevent copies
- [x] using Vec<u8> instead of `String` in the data structures
- [x] implement search
- [] implement streaming search
- [] when the minimizer is even, and its own reverse complement, how to break tie ? Using the superkmer ?
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
