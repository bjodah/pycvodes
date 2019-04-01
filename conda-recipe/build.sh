#!/bin/bash
export PYCVODES_LAPACK=openblas
# Sundials 2.7:
#export PYCVODES_SUNDIALS_LIBS=sundials_cvodes,sundials_nvecserial

# Sundials 3.1:
export PYCVODES_SUNDIALS_LIBS=sundials_cvodes,sundials_nvecserial,sundials_sunlinsollapackdense,sundials_sunlinsollapackband,sundials_sunlinsolklu,sundials_sunmatrixsparse

cat <<EOF>pycvodes/_config.py
env = {
    'LAPACK': "${PYCVODES_LAPACK}",
    'SUNDIALS_LIBS': "${PYCVODES_SUNDIALS_LIBS}",
    'NO_KLU': '0',
    'NO_LAPACK': '0',
    'SUNDIALS_PRECISION': 'double',
    'REAL_TYPE': 'double',
    'INDEX_TYPE': 'int32_t'
}
EOF
export PYCVODES_STRICT=1
export CPATH="${CPATH}:${PREFIX}/include"
python -m pip install --no-deps --ignore-installed .
