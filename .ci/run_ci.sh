#!/bin/bash -x
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -e

PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

export CPATH=${2}/include LIBRARY_PATH=${2}/lib LD_LIBRARY_PATH=${2}/lib  # SUNDIALS_ROOT=${2}
git clean -xfd

python3 setup.py sdist
(cd dist/; python3 -m pip install $PKG_NAME-$(python3 ../setup.py --version).tar.gz)
(cd /; python3 -m pytest --pyargs $PKG_NAME)
CXX=clang++-6.0 CC=clang-6.0 CFLAGS='-fsanitize=address' python3 -m pip install --force-reinstall .

PYTHONPATH=$(pwd) ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip
(cd /; python3 -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

cd tests/; make; make clean; cd -
cd tests/; make EXTRA_FLAGS=-DNDEBUG; make clean; cd -
cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -

(cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
(cd examples/; ../scripts/render_index.sh *.html)
