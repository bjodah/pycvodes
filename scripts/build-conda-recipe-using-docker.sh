#!/usr/bin/env bash
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
docker run -e HOST_USER_ID=$(id -u ${HOST_USER}) \
       -v "$REPOPATH":/mount \
       -v "$REPOPATH/.conda-cache/pkgs":/opt/conda/pkgs \
       -it condaforge/linux-anvil bash -c \
'set -xe; conda-build --output-folder /mount/dist /mount/conda-recipe'

#       -v "$REPOPATH/.conda-cache/conda-bld":/opt/conda/conda-bld \
