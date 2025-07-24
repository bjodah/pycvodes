#!/bin/bash
set -euxo pipefail
show_help() {
    echo "Usage: ./$(basename $0) [opts] /path/to/sundials/prefix/e.g./usr/local"
    echo "--native"
    echo "--python <python-exe-path>"
}

NATIVE=0
TEST_ASAN=0
MAKE_TMP_DIR=0
while [ $# -gt 1 ]; do
    case $1 in
        --native)
            NATIVE=1
            shift
            ;;
        --python)
            shift
            PYTHON=$1
            shift
            ;;
        --asan)
            TEST_ASAN=1
            shift
            ;;
        --tmp)
            MAKE_TMP_DIR=1
            shift
            ;;
        *)
            >&2 echo "Unrecognized parameter: $1"
            exit 1
    esac
done
if [ "$MAKE_TMP_DIR" = 1 ]; then
    REPO_TMP_DIR=$(mktemp -d)
    trap 'rm -rf -- "$REPO_TMP_DIR"' EXIT
    cp -ra . "$REPO_TMP_DIR/."
    cd "$REPO_TMP_DIR"
fi
SUNDBASE=$1
if [ ! -d "$SUNDBASE/include/sundials" ]; then >&2 echo "No such directory: $SUNDBASE"; exit 1; fi
if [[ $SUNDBASE =~ *-extended || $SUNDBASE =~ *-single ]]; then
    export PYCVODES_NO_LAPACK=1 PYCVODES_NO_KLU=1
fi
PYTHON=${PYTHON:-python3}
if [ $TEST_ASAN -eq 1 ]; then
    export PYTHON_ENV="env ASAN_OPTIONS=abort_on_error=1,detect_leaks=0"
else
    export PYTHON_ENV="env"
fi
LINKLIBS="$(${PYTHON_ENV} ${PYTHON} setup.py --print-linkline)"
export CPATH=/usr/include/suitesparse  # include <klu.h>
export CXXFLAGS="${CXXFLAGS:-} -isystem $SUNDBASE/include"
export LDFLAGS="$LINKLIBS -Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib -lopenblas"
export LD_LIBRARY_PATH=$(compgen -G "/opt-2/llvm-*/lib")

if [ $TEST_ASAN -eq 1 ]; then
    export CC=clang
    export CXX=clang++
    LLVM_ROOT=$(compgen -G "/opt-2/llvm-??")
    LLVM_LIB_DIR=$(compgen -G "${LLVM_ROOT}/lib/$(uname -m)-*")
    LIBCXX_ASAN_ROOT=$(compgen -G "/opt-2/libcxx*-asan")
    LIBCXX_ASAN_INCLUDE=${LIBCXX_ASAN_ROOT}/include/c++/v1
    if [ ! -d "$LIBCXX_ASAN_INCLUDE" ]; then
        >&2 echo "Not a directory: $LIBCXX_ASAN_INCLUDE"
        exit 1
    fi
    export CXXFLAGS="$CXXFLAGS -fsanitize=address -stdlib++-isystem ${LIBCXX_ASAN_INCLUDE} -ferror-limit=5"
    export LDFLAGS="${LDFLAGS:-} -fsanitize=address -Wl,-rpath,${LIBCXX_ASAN_ROOT}/lib -L${LIBCXX_ASAN_ROOT}/lib -lc++ -lc++abi -stdlib=libc++"
    export LIBRARY_PATH="$LLVM_ROOT/lib:${LIBCXX_ASAN_ROOT}/lib:${LIBRARY_PATH:-}"
    #LD_PRELOAD=

    export PYTHON_ENV="$PYTHON_ENV LD_PRELOAD=$(clang++ --print-file-name=libclang_rt.asan.so):${LIBCXX_ASAN_ROOT}/lib/libc++.so.1.0:${LIBCXX_ASAN_ROOT}/lib/libc++abi.so.1.0:${LIBCXX_ASAN_ROOT}/lib/libunwind.so.1.0"  # Or this failure appears:
    # AddressSanitizer: CHECK failed: asan_interceptors.cpp:463 "((__interception::real___cxa_throw)) != (0)" (0x0, 0x0)
    OPENMP_LIB="-Wl,-rpath,${LLVM_LIB_DIR} -lomp"
else
    export CC=gcc
    export CXX=g++
    export CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG -D_GLIBCXX_PEDANTIC"
    export CONTEXT="echo ''; echo ''; valgrind --error-exitcode=1"
    OPENMP_LIB="-lgomp"
fi


if [ -d ./build ]; then
    rm -r ./build
fi

CC=$CXX CFLAGS="$CXXFLAGS -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION" $PYTHON_ENV $PYTHON setup.py build_ext -i

export PYTHONPATH=$(pwd)

if [[ $SUNDBASE =~ .*-single ]]; then
    EXTRA_PYTEST_FLAGS="-k not test_get_include and not test_examples"
else
    EXTRA_PYTEST_FLAGS="-k not test_get_include"
fi

#gdb -ex r -args

# env \
#     LD_PRELOAD=\
# $(clang++ --print-file-name=libclang_rt.asan.so):\
# /opt-2/libcxx18-asan/lib/libc++.so.1:\
# /opt-2/libcxx18-asan/lib/libc++abi.so:\
# /opt-2/libcxx18-asan/lib/libunwind.so \
#gdb -ex r -args
$PYTHON_ENV $PYTHON -m pytest -sv "$EXTRA_PYTEST_FLAGS"


if [[ $SUNDBASE =~ .*-single ]]; then
    :
else
    cd tests/
    make clean
    make PYTHON="$PYTHON_ENV ${PYTHON}" OPENMP_LIB="${OPENMP_LIB}"
    cd -
fi
