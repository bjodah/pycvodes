#!/bin/bash -x
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -e

PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

for p in "${@:2}"
do
export CPATH=$p/include:$CPATH LIBRARY_PATH=$p/lib:$LIBRARY_PATH LD_LIBRARY_PATH=$p/lib:$LD_LIBRARY_PATH
done

git clean -xfd

python3 setup.py sdist
(cd dist/; python3 -m pip install $PKG_NAME-$(python3 ../setup.py --version).tar.gz)
(cd /; python3 -m pytest --pyargs $PKG_NAME)
(cd /; python3 -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

CXX=clang++-6.0 CC=clang-6.0 CFLAGS='-fsanitize=address' python3 -m pip install --force-reinstall .

if [[ "${LOW_PRECISION:-0}" != "1" ]]; then
    PYTHONPATH=$(pwd) ./scripts/run_tests.sh
    cd tests/; make; make clean; cd -
    cd tests/; make EXTRA_FLAGS=-DNDEBUG; make clean; cd -
    if [[ "${TEST_NATIVE_CLANG:-1}" == "1" ]]; then
        cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
        cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -
    fi

    (cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
    (cd examples/; ../scripts/render_index.sh *.html)
fi

