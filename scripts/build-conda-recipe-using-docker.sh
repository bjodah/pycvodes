#!/usr/bin/env bash
#
# Usage:
#
#  $ ./scripts/build-conda-recipe-using-docker.sh
# or
#  $ ./scripts/build-conda-recipe-using-docker.sh dist/conda-recipe-0.11.4
#
REPOPATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
if [ ! -d "$REPOPATH/dist" ]; then
    mkdir "$REPOPATH/dist"
fi
if [ ! -d "$REPOPATH/.conda-cache/pkgs" ]; then
    mkdir "$REPOPATH/.conda-cache/pkgs"
fi
# if [ ! -d "$REPOPATH/.conda-cache/conda-bld" ]; then
#     mkdir "$REPOPATH/.conda-cache/conda-bld"
# fi
HOST_USER=${SUDO_USER:-${LOGNAME}}
set -x
docker run \
       -e HOST_USER_ID=$(id -u ${HOST_USER}) \
       -v "$REPOPATH":/mount \
       -v "$REPOPATH/.conda-cache/pkgs":/opt/conda/pkgs \
       -it bjodah/anfilte bash -c \
       "set -xe; conda-build --output-folder /mount/dist /mount/${1:-conda-recipe} ${@:2}"
