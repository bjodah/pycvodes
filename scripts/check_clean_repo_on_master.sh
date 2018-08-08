#!/bin/bash
if [[ $(git rev-parse --abbrev-ref HEAD) != pycvodes-0.10.x ]]; then
    echo "We are not on the 'pycvodes-0.10.x' branch. Aborting..."
    exit 1
fi
if [[ ! -z $(git status -s) ]]; then
    echo "'git status' show there are some untracked/uncommited changes. Aborting..."
    exit 1
fi
