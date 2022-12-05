#!/bin/bash -x

#${PYTHON:-python3} -m pip uninstall -y pycvodes
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -eu
PKG_NAME=$1
SUNDBASE=${2%/}
set +u
if [ ! -e "$SUNDBASE/include/sundials/sundials_config.h" ]; then
    >&2 echo "Not a valid prefix for sundials: $SUNDBASE"
    exit 1
fi

export CFLAGS="-isystem $SUNDBASE/include $CFLAGS"
export LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib $LDFLAGS"

git clean -xfd

PKG_VERSION=$(${PYTHON:-python3} setup.py --version)
${PYTHON:-python3} setup.py sdist

cp -ra dist/ /tmp/
git clean -xfd


if [[ "${LOW_PRECISION:-0}" == "0" ]]; then
    if [ -d build/ ]; then rm -r build/; fi
    CXX=clang++-15 CC=clang-15 CFLAGS="-fsanitize=address -DPYCVODES_CLIP_TO_CONSTRAINTS=1 -UNDEBUG -O0 -g $CFLAGS" ${PYTHON:-python3} setup.py build_ext -i
    export PYTHON="env LD_PRELOAD=$(clang++-15 --print-file-name=libclang_rt.asan-$(uname -m).so) ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 ${PYTHON:-python3}"
    export PYTHONPATH=$(pwd)
    export ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-15/bin/llvm-symbolizer
    export ASAN_OPTIONS=symbolize=1
    gdb -ex r -args $PYTHON -m pytest -v
    
    LINKLIBS=$(${PYTHON} -c 'from pycvodes._libs import print_libs_linkline as pll; pll()')

    #unset LD_PRELOAD ASAN_OPTIONS

    cd tests/
    LDFLAGS="$LDFLAGS $LINKLIBS" make
    make CONTEXT="echo ''; echo ''; valgrind --error-exitcode=1"
    make clean

    LDFLAGS="$LDFLAGS $LINKLIBS" make EXTRA_FLAGS=-DNDEBUG
    make clean


    if [[ "${TEST_NATIVE_CLANG:-1}" == "1" ]]; then
        LDFLAGS="$LDFLAGS $LINKLIBS" LIBRARY_PATH=/usr/lib/llvm-15/lib:$LIBRARY_PATH make CXX=clang++-15 EXTRA_FLAGS=-fsanitize=address
        make clean
        
        LDFLAGS="$LDFLAGS $LINKLIBS" make CXX=clang++-15 EXTRA_FLAGS=-fsanitize=undefined
        make clean
    fi
    unset PYTHONPATH PYTHON
    cd -
    (cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
fi

if [[ "${BUILD_DOCS:-0}" == "1" ]]; then
    ${PYTHON:-python3} setup.py build_ext -i
    ${PYTHON:-python3} -m pytest --cov $PKG_NAME --cov-report html  # --pep8
    ./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
    ./scripts/generate_docs.sh
fi

(cd /tmp/dist/; CC=gcc-12 CXX=g++-12 ${PYTHON:-python3} -m pip install $PKG_NAME-$PKG_VERSION.tar.gz)
(cd /; ${PYTHON:-python3} -m pytest --verbose --pyargs $PKG_NAME)
(cd /; ${PYTHON:-python3} -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

