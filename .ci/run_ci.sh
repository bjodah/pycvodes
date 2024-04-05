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

${PYTHON:-python3} -m pip install cython pytest
PKG_VERSION=$(${PYTHON:-python3} setup.py --version)
${PYTHON:-python3} setup.py sdist

cp -ra dist/ /tmp/
git clean -xfd


if [[ "${LOW_PRECISION:-0}" == "0" ]]; then
    if [ -d build/ ]; then rm -r build/; fi
    # so one would think that the setuptools regression wrt. language='c++' would be solved 4 years
    # after the report (https://github.com/pypa/setuptools/issues/1732), but nope, as of 2024-04-03
    # it still isn't: https://github.com/pypa/distutils/pull/228 (but soon! ...maybe).
    #CXX=clang++ CC=clang
    LIBCXX_ASAN_ROOT=$(compgen -G "/opt-2/libcxx*-asan")
    env \
        CC=clang++ \
        CXX=clang++ \
        CFLAGS="-fsanitize=address -stdlib++-isystem ${LIBCXX_ASAN_ROOT}/include/c++/v1 -ferror-limit=5 -stdlib=libc++ -DPYCVODES_CLIP_TO_CONSTRAINTS=1 -UNDEBUG -O0 -g $CFLAGS" \
        LDFLAGS="-fsanitize=address -Wl,-rpath,${LIBCXX_ASAN_ROOT}/lib -L${LIBCXX_ASAN_ROOT}/lib -lc++ -lc++abi -stdlib=libc++ $LDFLAGS" \
        ${PYTHON:-python3} setup.py build_ext -i
    ORI_PYTHON="${PYTHON:-python3}"
    export PYTHON="env LD_PRELOAD=$(clang++ --print-file-name=libclang_rt.asan.so) ASAN_OPTIONS=abort_on_error=1,detect_leaks=0 ${PYTHON:-python3}"  # asan-$(uname -m).so
    ORI_PYTHONPATH="$PYTHONPATH"
    export PYTHONPATH=$(pwd)
    #export ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-15/bin/llvm-symbolizer
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
        LDFLAGS="$LDFLAGS $LINKLIBS" LIBRARY_PATH=$LLVM_ROOT/lib:$LIBRARY_PATH make CXX=clang++ EXTRA_FLAGS=-fsanitize=address
        make clean
        
        LDFLAGS="$LDFLAGS $LINKLIBS" make CXX=clang++ EXTRA_FLAGS=-fsanitize=undefined
        make clean
    fi
    export PYTHONPATH=$ORI_PYTHONPATH PYTHON=$ORI_PYTHON
    cd -
    (cd examples/; jupyter nbconvert --to=html --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=300 *.ipynb)
fi

if [[ "${BUILD_DOCS:-0}" == "1" ]]; then
    ${PYTHON:-python3} setup.py build_ext -i
    ${PYTHON:-python3} -m pytest --cov $PKG_NAME --cov-report html  # --pep8
    ./scripts/coverage_badge.py htmlcov/ htmlcov/coverage.svg
    ./scripts/generate_docs.sh
fi

(cd /tmp/dist/; CC=gcc-13 CXX=g++-13 ${PYTHON:-python3} -m pip install $PKG_NAME-$PKG_VERSION.tar.gz)
(cd /; ${PYTHON:-python3} -m pytest --verbose --pyargs $PKG_NAME)
(cd /; ${PYTHON:-python3} -c "from pycvodes import get_include as gi; import os; assert 'cvodes_cxx.pxd' in os.listdir(gi())")

