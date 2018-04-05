#!/bin/bash -xeu
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

python3 setup.py sdist
(cd dist/; python3 -m pip install $PKG_NAME-$(python3 ../setup.py --version).tar.gz)
(cd /; python3 -m pytest --pyargs $PKG_NAME)


PYTHONPATH=$(pwd) ./scripts/run_tests.sh --cov $PKG_NAME --cov-report html
./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg

# Make sure repo is pip installable from git-archive zip
git archive -o /tmp/$PKG_NAME.zip HEAD
python3 -m pip install --force-reinstall /tmp/$PKG_NAME.zip
(cd /; python3 -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

cd tests/; make EXTRA_LIBS="${EXTRA_LIBS:-}"; make clean; cd -
cd tests/; make EXTRA_LIBS="${EXTRA_LIBS:-}" EXTRA_FLAGS=-DNDEBUG; make clean; cd -
cd tests/; make EXTRA_LIBS="${EXTRA_LIBS:-}" CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
cd tests/; make EXTRA_LIBS="${EXTRA_LIBS:-}" CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -

(cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
(cd examples/; ../scripts/render_index.sh *.html)
