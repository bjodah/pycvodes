#!/bin/bash
CPLUS_INCLUDE_PATH=${PREFIX}/include ${PYTHON} setup.py build
${PYTHON} setup.py install
