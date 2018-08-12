# This file is replaced by setup.py in distributions for tagged releases
import logging
import pickle
import os
import warnings
import io
import sys
import tempfile
from contextlib import contextmanager


try:
    import appdirs
except ImportError:
    appdirs = None

pipes = None

if 'pytest' not in sys.modules:
    try:
        from wurlitzer import pipes
    except ImportError:
        pass

if sys.version_info[0] == 2:
    class TemporaryDirectory(object):
        def __init__(self):
            self.path = tempfile.mkdtemp()

        def __enter__(self):
            return self.path

        def __exit__(self, exc, value, tb):
            shutil.rmtree(self.path)
else:
    TemporaryDirectory = tempfile.TemporaryDirectory


def _warn(msg):
    if os.environ.get("PYCVODES_STRICT", '0') == '1':
        raise RuntimeError(msg)
    else:
        warnings.warn(msg)


def _compiles_ok(codestring):
    from distutils.ccompiler import new_compiler
    from distutils.sysconfig import customize_compiler
    from distutils.errors import CompileError
    ntf = tempfile.NamedTemporaryFile(suffix='.cpp', delete=False)
    ntf.write(codestring.encode('utf-8'))
    ntf.close()
    compiler = new_compiler()
    customize_compiler(compiler)
    out = ''
    try:
        if 'pytest' in sys.modules:
            compiler.compile([ntf.name])
        else:
            with capture_stdout_stderr() as out:
                compiler.compile([ntf.name])
    except CompileError:
        _ok = False
    except Exception:
        _ok = False
        _warn("Failed test compilation of '%s':\n %s" % (codestring, out))
    else:
        _ok = True

    os.unlink(ntf.name)
    return _ok, out


def _attempt_compilation():
    _math_ok, _math_out = _compiles_ok('#include <math.h>')
    if not _math_ok:
        _warn("Failed to include math.h: %s" % _math_out)

    _sundials_ok, _sundials_out = _compiles_ok('#include <sundials/sundials_config.h>')
    if not _sundials_ok:
        _warn("sundials not in include path, set e.g. $CPLUS_INCLUDE_PATH (%s):\n%s" %
              (os.environ.get('CPLUS_INCLUDE_PATH', ''), _sundials_out))

    _sun3_ok, _sun3_out = _compiles_ok("""
    #include <stdio.h>
    #include <sundials/sundials_config.h>
    #if SUNDIALS_VERSION_MAJOR >= 3
    #include <sunmatrix/sunmatrix_dense.h>
    #else
    #error "Sundials 2?"
    #endif
    """)

    _sun3, _lapack_ok = False, False
    if _sun3_ok:
        _lapack_ok, _lapack_out = _compiles_ok("""
    #include <stdio.h>
    #include <sundials/sundials_config.h>
    #if defined(SUNDIALS_BLAS_LAPACK)
    #  include <sunlinsol/sunlinsol_lapackband.h>
    #else
    #  error "Sundials 3+ was not configured to use lapack"
    #endif
    """)
        if not _lapack_ok:
            _warn("lapack not enabled in the sundials (>=3) distribtuion:\n%s" % _lapack_out)
        _sun3 = True
    else:
        _sun2_ok, _sun2_out = _compiles_ok("""
    #include <sundials/sundials_config.h>
    #if defined(SUNDIALS_PACKAGE_VERSION)   /* == 2.7.0 */
    #include <cvodes/cvodes_spgmr.h>
    #else
    #error "Unkown sundials version"
    #endif
    """)
        if _sun2_ok:
            _sun3 = False
            _lapack_ok, _lapack_out = _compiles_ok("""
    #include <sundials/sundials_config.h>
    #if defined(SUNDIALS_BLAS_LAPACK)
    #  include <cvodes/cvodes_lapack.h>
    #else
    #  error "Sundials 2 was not configured to use lapack"
    #endif
    """)
            if not _lapack_ok:
                _warn("lapack not enabled in the sundials (<3) distribution:\n%s" % _lapack_out)
        else:
            _warn("Unknown sundials version:\n%s" % _sun2_out)

    if 'PYCVODES_LAPACK' in os.environ:
        if os.environ['PYCVODES_LAPACK'] in ('', '0'):
            _lapack_ok = False

    return locals()

logger = logging.getLogger(__name__)

env = None
if appdirs:
    cfg = os.path.join(appdirs.user_config_dir('pycvodes-0.10'),
                       'python%d.%d-env.pkl' % sys.version_info[:2])
    if locals().get('_PYCVODES_IGNORE_CFG', 0) == 0:
        if os.path.exists(cfg) and os.path.getsize(cfg):
            with open(cfg, 'rb') as ifh:
                env = pickle.load(ifh)
        else:
            logger.info("Path: '%s' does not exist, will run test compilations" % cfg)
    else:
        logger.info("ignoring contents of '%s' (running from setup.py)" % cfg)
else:
    logger.info("appdirs not installed, will run test compilations")


if env is None:
    class Empty:
        pass

    _r = Empty()
    _r.__dict__.update(_attempt_compilation())

    env = {
        'LAPACK': 'blas,lapack' if _r._lapack_ok else '',
        'SUNDIALS_LIBS': 'sundials_cvodes,sundials_nvecserial' + (
            ',sundials_sunlinsollapackdense,sundials_sunlinsollapackband' if (_r._sun3 and _r._lapack_ok) else ((
                ',sundials_sunlinsoldense,sundials_sunlinsolband,sundials_sunlinsolspgmr,sundials_sunlinsolspbcgs,'
                'sundials_sunlinsolsptfqmr,sundials_sunmatrixdense,sundials_sunmatrixband'
            ) if (_r._sun3 and not _r._lapack_ok) else '')
        ),
    }
    if appdirs:
        cfg_dir = os.path.dirname(cfg)
        if not os.path.exists(cfg_dir):
            os.mkdir(cfg_dir)
        with open(cfg, 'wb') as ofh:
            pickle.dump(env, ofh)


for k, v in list(env.items()):
    env[k] = os.environ.get('%s_%s' % ('PYCVODES', k), v)
