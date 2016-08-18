#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Tested with Sundials 2.6.2

import io
import os
import shutil
import sys
from setuptools import setup
from setuptools.extension import Extension


pkg_name = 'pycvodes'

# Cythonize .pyx file if it exists (not in source distribution)
ext_modules = []
include_dirs = []


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)


USE_CYTHON = os.path.exists(_path_under_setup('pycvodes', '_cvodes_numpy.pyx'))

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import numpy as np
    LLAPACK = os.environ.get('LLAPACK', 'lapack')
    include_dirs = [np.get_include(), _path_under_setup('include')]
    ext = '.pyx' if USE_CYTHON else '.cpp'
    sources = ['pycvodes/_cvodes_numpy'+ext]
    ext_modules = [Extension('pycvodes._cvodes_numpy', sources)]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, include_path=['./include'])
    ext_modules[0].language = 'c++'
    ext_modules[0].extra_compile_args = ['-std=c++11']
    ext_modules[0].include_dirs = [_path_under_setup('include')]
    ext_modules[0].libraries += ['sundials_cvodes',
                                 os.environ.get('LLAPACK', 'lapack'),
                                 'sundials_nvecserial']


PYCVODES_RELEASE_VERSION = os.environ.get('PYCVODES_RELEASE_VERSION', '')

if os.environ.get('CONDA_BUILD', '0') == '1':
    # http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
    try:
        PYCVODES_RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = _path_under_setup(pkg_name, '_release.py')

if len(PYCVODES_RELEASE_VERSION) > 1:
    if PYCVODES_RELEASE_VERSION[0] != 'v':
        raise ValueError("PYCVODES_RELEASE_VERSION does not start with 'v'")
    TAGGED_RELEASE = True
    __version__ = PYCVODES_RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(open(release_py_path).read())

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

tests = [
    'pycvodes.tests',
]

with io.open(_path_under_setup(pkg_name, '__init__.py'), 'rt',
             encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
assert 10 < len(short_description) < 255
long_description = io.open(_path_under_setup('README.rst'),
                           encoding='utf-8').read()
assert len(long_description) > 100

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_description,
    classifiers=classifiers,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/bjodah/' + pkg_name,
    license='BSD',
    packages=[pkg_name] + tests,
    package_data={pkg_name: ['include/*.*']},
    install_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    setup_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    ext_modules=ext_modules,
    include_dirs=include_dirs
)

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set PYCVODES_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
