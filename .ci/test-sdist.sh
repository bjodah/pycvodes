#!/bin/bash
set -euxo pipefail
export PATH="/opt-2/gcc-13/bin:$PATH"

source /opt-3/cpython-v3.11-apt-deb/bin/activate
git clean -xfd
python3 setup.py sdist
cd ./dist/
CC=gcc pip install pycvodes-*.tar.gz
pip install pytest-flakes pytest-cov
./scripts/run_tests.sh --cov pycvodes --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
./scripts/generate_docs.sh
