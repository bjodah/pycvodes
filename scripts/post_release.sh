#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/post_release.sh v1.2.3 myserver githubuser
#
VERSION=${1#v}
SERVER=$2
GITHUBUSER=$3

./scripts/update-gh-pages.sh v$VERSION
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
SDIST_FILE=dist/${PKG}-$VERSION.tar.gz
if [[ ! -f "$SDIST_FILE" ]]; then
    >&2 echo "Nonexistent file $SDIST_FILE"
    exit 1
fi
SHA256=$(openssl sha256 "$SDIST_FILE" | cut -f2 -d' ')
if [[ -d "dist/conda-recipe-$VERSION" ]]; then
    rm -r "dist/conda-recipe-$VERSION"
fi
cp -r conda-recipe/ dist/conda-recipe-$VERSION
sed -i -E \
    -e "s/\{\% set version(.+)/\{\% set version = \"$VERSION\" \%\}\n\{\% set sha256 = \"$SHA256\" \%\}/" \
    -e "s/git_url:(.+)/fn: \{\{ name \}\}-\{\{ version \}\}.tar.gz\n  url: https:\/\/pypi.io\/packages\/source\/\{\{ name\[0\] \}\}\/\{\{ name \}\}\/\{\{ name \}\}-\{\{ version \}\}.tar.gz\n  sha256: \{\{ sha256 \}\}/" \
    -e "/cython/d" \
    dist/conda-recipe-$VERSION/meta.yaml

for CONDA_PY in 3.5 3.6; do
    ./scripts/build-conda-recipe-using-docker.sh dist/conda-recipe-$VERSION  --python ${CONDA_PY}
done

scp dist/${PKG}*${VERSION}*.bz2 $PKG@$SERVER:~/public_html/conda-packages/
scp -r dist/conda-recipe-$VERSION/ $PKG@$SERVER:~/public_html/conda-recipes/
scp "$SDIST_FILE" "$PKG@$SERVER:~/public_html/releases/"
