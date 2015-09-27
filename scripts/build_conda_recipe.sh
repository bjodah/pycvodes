#!/bin/bash -uex
# Usage:
#
#    $ ./scripts/build_conda_recipe.sh v1.2.3
#
echo $1>__conda_version__.txt
trap "rm __conda_version__.txt" EXIT SIGINT SIGTERM
conda build conda-recipe
