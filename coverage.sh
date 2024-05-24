# ---------------------------------------------------------------------------- #
#
# Generate test coverage:
# * report on terminal
# * lcov file
# * html
#
# See CONTRIBUTING.md for dependency installation
#
# USAGE
# -----
# ./coverage.sh
# or
# ./coverage.sh <name_of_the_module_to_test>
#
# ---------------------------------------------------------------------------- #
module_name=$1

LCOV_FILE=.lcov.info
#
# Clean files
#
rm $LCOV_FILE 2>/dev/null
cargo llvm-cov clean --workspace
#
# Terminal visual report
#
cargo llvm-cov nextest $module_name
#
# Generate lcov file in .lcov.info
#
cargo llvm-cov report --lcov --output-path $LCOV_FILE
#
# Generate html file in target/llvm-cov/html
#
cargo llvm-cov report --html
