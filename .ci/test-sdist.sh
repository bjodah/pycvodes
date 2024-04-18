#!/bin/bash
set -euxo pipefail
export PATH="/opt-2/gcc-13/bin:$PATH"

source /opt-3/cpython-v3.11-apt-deb/bin/activate
git clean -xfd
python3 setup.py sdist
cd ./dist/
SUNDBASE=/opt-3/sundials-6.7.0-release
env \
    CFLAGS="-isystem $SUNDBASE/include -isystem /usr/include/suitesparse" \
    LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib" \
    CC=gcc CXX=g++ pip install ${CI_REPO_NAME}-*.tar.gz
pip install pytest-flakes pytest-cov matplotlib sphinx numpydoc sphinx-rtd-theme
../scripts/run_tests.sh --cov ${CI_REPO_NAME} --cov-report html
../scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

cd -
pip install .
./scripts/generate_docs.sh
