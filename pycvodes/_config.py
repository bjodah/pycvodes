# This file is replaced by setup.py in distributions for tagged releases
import os
import warnings
import io
import sys
import tempfile
from contextlib import contextmanager

pipes = None

if not 'pytest' in sys.modules:
    try:
        from wurlitzer import pipes
    except ImportError:
        pass

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
        if pipes is None:
            compiler.compile([ntf.name])
        else:
            with pipes() as out_err:
                compiler.compile([ntf.name])
            out = '\n'.join([p.read() for p in out_err])
    except CompileError:
        _ok = False
    except Exception:
        _ok = False
        _warn("Failed test compilation of '%s':\n %s" % (codestring, out))
    else:
        _ok = True

    os.unlink(ntf.name)
    return _ok, out


_math_ok, _math_out = _compiles_ok('#include <math.h>')
if not _math_ok:
    _warn("Failed to include math.h: %s" % _math_out)
    _warn(_math_out)

_sundials_ok, _sundials_out = _compiles_ok('#include <sundials/sundials_config.h>')
if not _sundials_ok:
    _warn("sundials not in include path, set e.g. $CPLUS_INCLUDE_PATH (%s):\n%s" %
          (os.environ.get('CPLUS_INCLUDE_PATH', ''), _sundials_out))

_sun3_ok, _sun3_out = _compiles_ok("""
#include <sundials/sundials_config.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <stdio.h>
#include <sunmatrix/sunmatrix_dense.h>
#else
#error "Sundials 2?"
#endif
""")

_sun3, _lapack_ok = False, False
if _sun3_ok:
    _lapack_ok, _lapack_out = _compiles_ok("""
#include <stdio.h>
#include <sunlinsol/sunlinsol_lapackband.h>""")
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
#include <cvodes/cvodes_lapack.h>

""")
        if not _lapack_ok:
            _warn("lapack not enabled in the sundials (<3) distribution:\n%s" % _lapack_out)
    else:
        _warn("Unknown sundials version:\n%s" % _sun2_out)

if 'PYCVODES_LAPACK' in os.environ:
    if os.environ['PYCVODES_LAPACK'] in ('', '0'):
        _lapack_ok = False

env = {
    'LAPACK': 'blas,lapack' if _lapack_ok else '',
    'SUNDIALS_LIBS': 'sundials_cvodes,sundials_nvecserial' + (
        ',sundials_sunlinsollapackdense,sundials_sunlinsollapackband' if (_sun3 and _lapack_ok) else ((
            ',sundials_sunlinsoldense,sundials_sunlinsolband,sundials_sunlinsolspgmr,sundials_sunlinsolspbcgs,'
            'sundials_sunlinsolsptfqmr,sundials_sunmatrixdense,sundials_sunmatrixband'
        ) if (_sun3 and not _lapack_ok) else '')
    ),
}

for k, v in list(env.items()):
    env[k] = os.environ.get('%s_%s' % ('PYCVODES', k), v)
