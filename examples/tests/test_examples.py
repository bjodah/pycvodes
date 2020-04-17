# -*- coding: utf-8 -*-

import glob
import os
import subprocess
import sys

import pytest


tests = glob.glob(os.path.join(os.path.dirname(__file__), '../*.py'))


@pytest.mark.parametrize('pypath', tests)
def test_examples(pypath):
    p = subprocess.Popen([sys.executable, pypath])
    assert p.wait() == 0  # SUCCESS==0

    p = subprocess.Popen([sys.executable, pypath, '--nt', '2'])
    assert p.wait() == 0  # SUCCESS==0
