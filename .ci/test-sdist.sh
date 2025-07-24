#!/bin/bash
set -euxo pipefail
#export PATH="/opt-2/gcc-13/bin:$PATH"

source /opt-3/cpython-v3.*-apt-deb/bin/activate
git clean -xfd
python3 setup.py sdist
cd ./dist/
SUNDBASE=$(compgen -G "/opt-3/sundials-6.*-release")
export CFLAGS="-isystem $SUNDBASE/include -isystem /usr/include/suitesparse"
env \
    CXXFLAGS="${CFLAGS}" \
    LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib" \
    CC=gcc CXX=g++ pip install ${CI_REPO_NAME}-*.tar.gz
python3 -m pytest --pyargs pycvodes
cd -
env \
    CXXFLAGS="${CFLAGS}" \
    LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib" \
    CC=gcc CXX=g++ python3 setup.py build_ext -i

pip install pytest-flakes pytest-cov matplotlib sphinx numpydoc sphinx-rtd-theme
./scripts/run_tests.sh --cov ${CI_REPO_NAME} --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
./scripts/generate_docs.sh
