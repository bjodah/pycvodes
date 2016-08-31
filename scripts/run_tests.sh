#!/bin/bash -e
# Usage
#   $ ./scripts/run_tests.sh
# or
#   $ ./scripts/run_tests.sh --cov pycvodes --cov-report html
ulimit -v 524288  # 512 MiB virt. mem: needed for testing MemoryError (std::bad_alloc)
${PYTHON:-python} setup.py build_ext -i
${PYTHON:-python} -m pytest --doctest-modules --pep8 --flakes $@
${PYTHON:-python} -m doctest README.rst
