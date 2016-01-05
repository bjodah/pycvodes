#!/bin/bash
touch doc/_build/html/.nojekyll
cp LICENSE doc/_build/html/.nojekyll
git config --global user.name "drone"
git config --global user.email "drone@nohost.com"
mkdir -p deploy/public_html/branches/"${CI_BRANCH}" deploy/script_queue
cp -r dist/* htmlcov/ examples/ doc/_build/html/ deploy/public_html/branches/"${CI_BRANCH}"/
if bash -c '[[ "$CI_BRANCH" == "master" ]]'; then
    mkdir -p deploy/gh-pages
    ./scripts/dir_to_branch.sh doc/_build/html bjodah "${CI_REPO##*/}" gh-pages deploy/gh-pages
fi
