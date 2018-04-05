#!/bin/bash -eu
PREFIX=$1
if [ ! -d "$PREFIX" ]; then 2>&1 echo "Not a directory: $PREFIX"; exit 1; fi
wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh -O ~/miniconda.sh
/bin/bash ~/miniconda.sh -b -p "$PREFIX"
rm ~/miniconda.sh
PATH="$PREFIX/bin:$PATH" conda install conda-build
