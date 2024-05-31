#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
${PYTHON:-python3} -m pytest --pyargs pycvodes --doctest-modules --flakes "$@"
MPLBACKEND=Agg ${PYTHON:-python3} -m doctest $(dirname $(realpath $0))/../README.rst
