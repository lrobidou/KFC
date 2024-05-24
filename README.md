# KFC

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
