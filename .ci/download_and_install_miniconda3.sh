#!/bin/bash -eu
PREFIX=$1
if [ -d "$PREFIX" ]; then 2>&1 echo "Directory already exists: $PREFIX"; exit 1; fi
wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh -O ~/miniconda.sh
/bin/bash ~/miniconda.sh -b -p "$PREFIX"
rm ~/miniconda.sh
PATH="$PREFIX/bin:$PATH" conda config --set always_yes yes && \
PATH="$PREFIX/bin:$PATH" conda config --set changeps1 no && \
PATH="$PREFIX/bin:$PATH" conda config --set show_channel_urls True && \
PATH="$PREFIX/bin:$PATH" conda config --add channels conda-forge && \
PATH="$PREFIX/bin:$PATH" conda install --quiet conda-build
