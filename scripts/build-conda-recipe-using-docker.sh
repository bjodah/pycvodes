#!/usr/bin/env bash
#
# Usage:
#
#  $ ./scripts/build-conda-recipe-using-docker.sh conda/recipe
# or
#  $ ./scripts/build-conda-recipe-using-docker.sh dist/conda-recipe-0.11.4
#
REPOPATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
if [ ! -d "$REPOPATH/dist" ]; then
    mkdir "$REPOPATH/dist"
fi
CMD="set -xe; conda-build "${@:2}" --output-folder /mount/dist /mount/$1"
set -x
CONTAINERNAME=anfilte-$(basename $REPOPATH)
HOST_USER=${SUDO_USER:-${LOGNAME}}

# TODO: using docker commit it should be possible to use containers as cache
docker run --name $CONTAINERNAME \
       -e HOST_USER_ID=$(id -u ${HOST_USER}) \
       -v "$REPOPATH":/mount \
       -it bjodah/anfilte \
       bash -c "$CMD"

EXIT_CODE=$(docker inspect ${CONTAINERNAME} --format='{{.State.ExitCode}}')
docker rm $CONTAINERNAME
exit ${EXIT_CODE}
