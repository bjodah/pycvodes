#!/bin/bash -x
set -euxo pipefail
export PATH=/opt-2/gcc-13/bin:$PATH
export CC=gcc
export CXX=g++
#                                     .ci/test-case.sh --asan --python /opt-3/cpython-v3.11.9-asan/bin/python3 /opt-3/sundials-6.7.0-asan
                                       .ci/test-case.sh        --python /opt-3/pypy3.10-v7.3.15-linux64/bin/python3 /opt-3/sundials-6.7.0-debug
                                       .ci/test-case.sh        --python /opt-3/cpython-v3.11.?-debug/bin/python3 /opt-3/sundials-6.7.0-debug
                                       .ci/test-case.sh        --python /opt-3/cpython-v3.12.2-debug/bin/python3 /opt-3/sundials-6.7.0-debug
PYCVODES_NO_LAPACK=1 PYCVODES_NO_KLU=1 .ci/test-case.sh        --python /opt-3/cpython-v3.11-apt-deb/bin/python3 /opt-3/sundials-6.7.0-extended
PYCVODES_NO_LAPACK=1 PYCVODES_NO_KLU=1 .ci/test-case.sh        --python /opt-3/cpython-v3.12.2-release/bin/python3 /opt-3/sundials-6.7.0-single
