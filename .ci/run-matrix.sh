#!/bin/bash -x
set -euxo pipefail

.ci/test-case.sh --asan --python /opt-3/cpython-v3.11.9-asan/bin/python3 /opt-3/sundials-6.7.0-asan
.ci/test-case.sh        --python /opt-3/cpython-v3.12.2-debug/bin/python3 /opt-3/sundials-6.7.0-debug
.ci/test-case.sh        --python /opt-3/cpython-v3.11-apt-deb/bin/python3 /opt-3/sundials-6.7.0-extended
.ci/test-case.sh        --python /opt-3/cpython-v3.12.2-release/bin/python3 /opt-3/sundials-6.7.0-single
