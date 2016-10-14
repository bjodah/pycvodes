#!/bin/bash
LLAPACK=openblas CPLUS_INCLUDE_PATH=${PREFIX}/include ${PYTHON} setup.py build
${PYTHON} setup.py install --single-version-externally-managed --record record.txt
