#!/bin/bash -xeu
# Usage:
#
#    $ ./scripts/post_release.sh v1.2.3 myserver githubuser
#
VERSION=${1#v}
SERVER=$2
GITHUBUSER=$3
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
    -e "s/\{\% set version(.+)/\{\% set version = \"$VERSION\" \%\}/" \
    -e "s/git_url:(.+)/fn: $PKG-$VERSION.tar.gz\n  url: https:\/\/pypi.io\/packages\/source\/${PKG:0:1}\/${PKG}\/${PKG}-$VERSION.tar.gz\n  sha256: $SHA256/" \
    -e "/cython/d" \
    dist/conda-recipe-$VERSION/meta.yaml

env ${PKG_UPPER}_RELEASE_VERSION=v$VERSION python setup.py upload_sphinx

# Specific for this project:
scp -r dist/conda-recipe-$VERSION/ $PKG@$SERVER:~/public_html/conda-recipes/
scp "$SDIST_FILE" "$PKG@$SERVER:~/public_html/releases/"
for CONDA_PY in 2.7 3.4 3.5; do
    for CONDA_NPY in 1.11; do
        ssh $PKG@$SERVER "source /etc/profile; conda-build --python $CONDA_PY --numpy $CONDA_NPY ~/public_html/conda-recipes/conda-recipe-$VERSION/"
    done
done
