#!/bin/bash -xe
# Usage:
#
#    $ ./scripts/release.sh v1.2.3 ~/anaconda2/bin
#

if [[ $1 != v* ]]; then
    echo "Argument does not start with 'v'"
    exit 1
fi
./scripts/check_clean_repo_on_master.sh
cd $(dirname $0)/..
# PKG will be name of the directory one level up containing "__init__.py" 
PKG=$(find . -maxdepth 2 -name __init__.py -print0 | xargs -0 -n1 dirname | xargs basename)
PKG_UPPER=$(echo $PKG | tr '[:lower:]' '[:upper:]')
./scripts/run_tests.sh
env ${PKG_UPPER}_RELEASE_VERSION=$1 python setup.py sdist
PATH=$2:$PATH ./scripts/build_conda_recipe.sh $1 27 34

# All went well, add a tag and push it.
git tag -a $1 -m $1
git push
git push --tags
VERSION=${1#v}
twine upload dist/${PKG}-$VERSION.tar.gz
MD5=$(md5sum dist/${PKG}-$VERSION.tar.gz | cut -f1 -d' ')

if [[ -d dist/conda-recipe-${1#v} ]]; then
    rm -r dist/conda-recipe-${1#v}
fi
cp -r conda-recipe/ dist/conda-recipe-${1#v}
sed -i -E -e "s/version:(.+)/version: $VERSION/" -e "s/path:(.+)/fn: $PKG-$VERSION.tar.gz\n    url: https:\/\/pypi.python.org\/packages\/source\/${PKG:0:1}\/$PKG\/$PKG-$VERSION.tar.gz#md5=$MD5\n    md5: $MD5/" dist/conda-recipe-${1#v}/meta.yaml

# Specific for this project:
SERVER=hera
scp -r dist/conda-recipe-${1#v}/ $PKG@$SERVER:~/public_html/conda-recipes/
scp dist/${PKG}-$VERSION.tar.gz $PKG@$SERVER:~/public_html/releases/
ssh $PKG@$SERVER "source /etc/profile; CONDA_PY=27 conda-build ~/public_html/conda-recipes/conda-recipe-${1#v}/"
ssh $PKG@$SERVER "source /etc/profile; CONDA_PY=34 conda-build ~/public_html/conda-recipes/conda-recipe-${1#v}/"

