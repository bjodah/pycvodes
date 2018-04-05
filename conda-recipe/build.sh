#!/bin/bash
sed -i -E -e 's/'"'"'lapack'"'"'/'"'"'openblas'"'"'/' pycvodes/_config.py
export PYCVODES_STRICT=1
export PYCVODES_LAPACK=openblas
export CPATH=${PREFIX}/include
export LIBRARY_PATH=${PREFIX}/lib
export LD_LIBRARY_PATH=${PREFIX}/lib
${PYTHON} -m pip install --no-deps --ignore-installed .
