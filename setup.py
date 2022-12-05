#!/usr/bin/env python
# -*- coding: utf-8 -*-


import io
import logging
import os
import re
import shutil
import subprocess
import sys
import warnings
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
try:
    import cython
except ImportError:
    _HAVE_CYTHON = False
else:
    _HAVE_CYTHON = True
    assert cython  # silence pep8


pkg_name = 'pycvodes'
url = 'https://github.com/bjodah/' + pkg_name
license = 'BSD'


def _path_under_setup(*args):
    return os.path.join(*args)

release_py_path = _path_under_setup(pkg_name, '_release.py')


_version_env_var = '%s_RELEASE_VERSION' % pkg_name.upper()
RELEASE_VERSION = os.environ.get(_version_env_var, '')

if len(RELEASE_VERSION) > 1:
    if RELEASE_VERSION[0] != 'v':
        raise ValueError("$%s does not start with 'v'" % _version_env_var)
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
else:  # set `__version__` from _release.py:
    TAGGED_RELEASE = False
    with open(release_py_path) as ifh:
        exec(ifh.read())
    if __version__.endswith('git'):
        try:
            _git_version = subprocess.check_output(
                ['git', 'describe', '--dirty']).rstrip().decode('utf-8')
        except subprocess.CalledProcessError:
            warnings.warn("A git-archive is being installed - version information incomplete.")
        else:
            if 'develop' not in sys.argv:
                warnings.warn("Using git to derive version: dev-branches may compete.")
                _ver_tmplt = r'\1.post\2' if os.environ.get('CONDA_BUILD', '0') == '1' else r'\1.post\2+\3'
                __version__ = re.sub(r'v([0-9.]+)-(\d+)-(\S+)', _ver_tmplt, _git_version)  # .dev < '' < .post

package_include = os.path.join(pkg_name, 'include')

_src = {ext: _path_under_setup(pkg_name, '_cvodes.' + ext) for ext in "cpp pyx".split()}
if _HAVE_CYTHON and os.path.exists(_src["pyx"]):
    # Possible that a new release of Python needs a re-rendered Cython source,
    # or that we want to include possible bug-fix to Cython, disable by manually
    # deleting .pyx file from source distribution.
    USE_CYTHON = True
    if os.path.exists(_src['cpp']):
        os.unlink(_src['cpp'])  # ensure c++ source is re-generated.
else:
    USE_CYTHON = False

ext_modules = []

if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import numpy as np
    env = None  # silence pyflakes, 'env' is actually set on the next line
    _PYCVODES_IGNORE_CFG = 1  # avoid using cached config upon running setup.py
    env = dict(
        LAPACK="blas,lapack",
        NO_KLU="0",
        NO_LAPACK="0",
        SUNDIALS_LIBS="",
    )
    for k, v in list(env.items()):
        env[k] = os.environ.get('%s_%s' % (pkg_name.upper(), k), v)
    _USE_KLU = env.get('NO_KLU', '0') == '0'
    if env.get('NO_LAPACK', '0') == '1' or env['LAPACK'] in ('', '0'):
        _USE_LAPACK = False
    else:
        _USE_LAPACK = True

    if env["SUNDIALS_LIBS"] == "":
        get_libs = None  # silence linting check (pyflakes), get from _libs.py
        exec(open(_path_under_setup(pkg_name, "_libs.py")).read())
        env["SUNDIALS_LIBS"] = get_libs(dict(LAPACK=_USE_LAPACK, KLU=_USE_KLU))

    logger = logging.getLogger(__name__)
    logger.info("Config for pycvodes: %s" % str(env))
    sources = [_src["pyx" if USE_CYTHON else "cpp"]]
    ext_modules = [Extension('%s._cvodes' % pkg_name, sources)]
    if USE_CYTHON:
        from Cython.Build import cythonize
        ext_modules = cythonize(ext_modules, include_path=[
            package_include,
            os.path.join('external', 'anyode', 'cython_def')
        ])
    ext_modules[0].language = 'c++'
    ext_modules[0].include_dirs = [np.get_include(), package_include,
                                   os.path.join('external', 'anyode', 'include')]

    ext_modules[0].define_macros += [
        ('PYCVODES_NO_KLU', '0' if _USE_KLU else '1'),
        ('PYCVODES_NO_LAPACK', '0' if _USE_LAPACK else '1'),
        ('ANYODE_NO_LAPACK', '0' if _USE_LAPACK else '1')
    ]

    if env['SUNDIALS_LIBS']:
        ext_modules[0].libraries += env['SUNDIALS_LIBS'].split(',')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append('--std=c++14')
            if sys.platform == 'darwin' and re.search("clang", self.compiler.compiler[0]) is not None:
                opts += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


classifiers = [
    "Development Status :: 4 - Beta",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

tests = [
    '%s.tests' % pkg_name,
]

with io.open(_path_under_setup(pkg_name, '__init__.py'), 'rt', encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
if not 10 < len(short_description) < 255:
    warnings.warn("Short description from __init__.py proably not read correctly.")
long_description = io.open(_path_under_setup('README.rst'),
                           encoding='utf-8').read()
if not len(long_description) > 100:
    warnings.warn("Long description from README.rst probably not read correctly.")
_author, _author_email = io.open(_path_under_setup('AUTHORS'), 'rt', encoding='utf-8').readline().split('<')

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_description,
    classifiers=classifiers,
    author=_author.strip(),
    author_email=_author_email.split('>')[0].strip(),
    url=url,
    license=license,
    packages=[pkg_name] + tests,
    include_package_data=True,
    install_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    setup_requires=['numpy'] + (['cython'] if USE_CYTHON else []),
    extras_require={'docs': ['Sphinx', 'sphinx_rtd_theme', 'numpydoc']},
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt}
)

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist files
            # depending on tagged version (set PYCVODES_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            with open(release_py_path, 'wt') as ofh:
                ofh.write("__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
