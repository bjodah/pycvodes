#!/bin/bash -x

#${PYTHON:-python3} -m pip uninstall -y pycvodes
rm -r /usr/local/lib/python*/dist-packages/pycvodes*  # pip uninstall is useless
set -eu
PKG_NAME=$1
SUNDBASE=$2
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
(cd /tmp/dist/; CC=gcc-11 CXX=g++-11 ${PYTHON:-python3} -m pip install $PKG_NAME-$PKG_VERSION.tar.gz)
(cd /; ${PYTHON:-python3} -m pytest --verbose --pyargs $PKG_NAME)
(cd /; ${PYTHON:-python3} -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

if [[ "${LOW_PRECISION:-0}" != "1" ]]; then
    (cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
    (cd examples/; ../scripts/render_index.sh *.html)
fi

if [[ "${BUILD_DOCS:-0}" == "1" ]]; then
    ${PYTHON:-python3} setup.py build_ext -i
    ${PYTHON:-python3} -m pytest --cov $PKG_NAME --cov-report html  # --pep8
    ./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
    ${PYTHON:-python3} setup.py sdist
    ./scripts/generate_docs.sh
fi

if [[ "${LOW_PRECISION:-0}" != "1" ]]; then
    if [ -d build/ ]; then rm -r build/; fi
    CXX=clang++-12 CC=clang-12 CFLAGS="-fsanitize=address -DPYCVODES_CLIP_TO_CONSTRAINTS=1 $CFLAGS" ${PYTHON:-python3} setup.py build_ext -i
    export PYTHON="env LD_PRELOAD=$(clang++-12 -print-file-name libclang_rt.asan-$(uname -m).so) ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 ${PYTHON:-python3}" PYTHONPATH=$(pwd)
    export ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-12/bin/llvm-symbolizer
    export ASAN_OPTIONS=symbolize=1
    ./scripts/run_tests.sh
    LINKLIBS=$(${PYTHON} -c "from pycvodes._libs import print_libs_linkline as pll; pll()")
    #unset LD_PRELOAD ASAN_OPTIONS
    cd tests/; LDFLAGS="$LDFLAGS $LINKLIBS" make; make clean; cd -
    cd tests/; LDFLAGS="$LDFLAGS $LINKLIBS" make EXTRA_FLAGS=-DNDEBUG; make clean; cd -
    if [[ "${TEST_NATIVE_CLANG:-1}" == "1" ]]; then
        cd tests/; LDFLAGS="$LDFLAGS $LINKLIBS" LIBRARY_PATH=/usr/lib/llvm-12/lib:$LIBRARY_PATH make CXX=clang++-12 EXTRA_FLAGS=-fsanitize=address; make clean; cd -
        cd tests/; LDFLAGS="$LDFLAGS $LINKLIBS" make CXX=clang++-12 EXTRA_FLAGS=-fsanitize=undefined; make clean; cd -
    fi
    unset PYTHONPATH PYTHON
fi
