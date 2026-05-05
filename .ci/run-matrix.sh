#!/bin/bash -x
llvm_env_sh=$(compgen -G "/etc/profile.d/llvm-??.sh")
set -euxo pipefail
export CC=gcc
export CXX=g++
if [ -e $llvm_env_sh ]; then source $llvm_env_sh ; fi
                                       .ci/test-case.sh --asan --python /opt-3/cpython-v3.13*-asan/bin/python3      /opt-3/sundials-6.7.0-asan
                                       .ci/test-case.sh        --python /opt-3/pypy3.1*-linux64/bin/python3         /opt-3/sundials-6.7.0-debug
                                       .ci/test-case.sh        --python /opt-3/cpython-v3.11.*-release/bin/python3  /opt-3/sundials-6.7.0-debug
                                       .ci/test-case.sh        --python /opt-3/cpython-v3.12.*-debug/bin/python3    /opt-3/sundials-6.7.0-debug
PYCVODES_NO_LAPACK=1 PYCVODES_NO_KLU=1 .ci/test-case.sh        --python /opt-3/cpython-v3.13.*-release/bin/python3  /opt-3/sundials-6.7.0-single

source /opt-3/cpython-v3.*-apt-deb/bin/activate
PYCVODES_NO_LAPACK=1 PYCVODES_NO_KLU=1 .ci/test-case.sh        --python python                                      /opt-3/sundials-6.7.0-extended
