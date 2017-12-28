#!/bin/bash -xeu
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

python2.7 setup.py sdist
for PYTHON in python2.7 python3; do
    (cd dist/; $PYTHON -m pip install $PKG_NAME-$($PYTHON ../setup.py --version).tar.gz)
    (cd /; $PYTHON -m pytest --pyargs $PKG_NAME)
done

PYTHONPATH=$(pwd) ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip
(cd /; python3 -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

cd tests/; make; cd -
cd tests/; make EXTRA_FLAGS=-DNDEBUG; cd -
cd tests/; make CXX=clang++-5.0 EXTRA_FLAGS=-fsanitize=address -B; cd -
cd tests/; make CXX=clang++-5.0 EXTRA_FLAGS=-fsanitize=undefined -B; cd -
