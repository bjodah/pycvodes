#!/bin/bash -x
set -euxo pipefail

for flags in "" "--asan"; do
    .ci/test-config.sh $flags --native /opt-3/sundials-6.7.0-release
    .ci/test-config.sh $flags --native /opt-3/sundials-6.7.0-debug
    .ci/test-config.sh $flags --native /opt-3/sundials-6.7.0-single
    .ci/test-config.sh $flags --native /opt-3/sundials-6.7.0-extended
    .ci/test-config.sh $flags --python /opt-3/cpython-v3.11.9-debug/bin/python3 /opt-3/sundials-6.7.0-debug
    .ci/test-config.sh $flags --python /opt-3/cpython-v3.12.2-debug/bin/python3 /opt-3/sundials-6.7.0-single
    .ci/test-config.sh $flags --python /opt-3/cpython-v3.11-apt-deb/bin/python3 /opt-3/sundials-6.7.0-extended
    .ci/test-config.sh $flags --python /opt-3/cpython-v3.12.2-release/bin/python3 /opt-3/sundials-6.7.0-release
done

