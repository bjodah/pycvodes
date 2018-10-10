#!/bin/bash -x

set +e
#python3 -m pip uninstall -y pycvodes
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -e

PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${PKG_NAME^^}_RELEASE_VERSION=\$CI_BRANCH
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

set +e
#python3 -m pip uninstall -y pycvodes
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -e

if [ -d build/ ]; then rm -r build/; fi
CXX=clang++-6.0 CC=clang-6.0 CFLAGS='-fsanitize=address' python3 setup.py build_ext -i

if [[ "${LOW_PRECISION:-0}" != "1" ]]; then
    PYTHONPATH=$(pwd) ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 LD_PRELOAD=/usr/lib/llvm-6.0/lib/clang/6.0.1/lib/linux/libclang_rt.asan-x86_64.so ./scripts/run_tests.sh
    cd tests/; make; make clean; cd -
    cd tests/; make EXTRA_FLAGS=-DNDEBUG; make clean; cd -
    if [[ "${TEST_NATIVE_CLANG:-1}" == "1" ]]; then
        cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
        cd tests/; make CXX=clang++-6.0 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -
    fi

    (cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
    (cd examples/; ../scripts/render_index.sh *.html)
fi

