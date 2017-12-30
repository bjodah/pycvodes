import os
import subprocess


def _compiles_ok(codestring):
    _preproc = subprocess.Popen('cc -E -x c++ -', shell=True, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _err = _preproc.communicate(codestring.encode('utf-8'))
    _retcode = _preproc.wait()
    if _retcode == 0:
        return True, out
    elif _retcode > 0:
        return False, out
    else:
        raise RuntimeError("Failed test compilation of '%s':\n%s" % (codestring, out))

_dense_ok, _dense_out = _compiles_ok('#include <sunmatrix/sunmatrix_dense.h>')
if _dense_ok:
    _sun3 = True
    _lapack_ok, _lapack_out = _compiles_ok('#include <sunlinsol/sunlinsol_lapackband.h>')
    if not _lapack_ok:
        raise RuntimeError("lapack not enabled in the sundials (>=3) distribtuion:\n%s" % _lapack_out)
else:
    _spgmr_ok, _spgmr_out = _compiles_ok('#include <cvodes/cvodes_spgmr.h>')
    if _spgmr_ok:
        _sun3 = False
        _lapack_ok, _lapack_out = _compiles_ok('#include <cvodes/cvodes_lapack.h>')
        if not _lapack_ok:
            raise RuntimeError("lapack not enabled in the sundials (<3) distribution:\n%s" % _lapack_out)
    else:
        raise RuntimeError("sundials not in include path, set e.g. $CPLUS_INCLUDE_PATH (%s):\n%s" %
                           (os.environ.get('CPLUS_INCLUDE_PATH', ''), _spgmr_out))

env = {
    'LAPACK': 'lapack',
    'SUNDIALS_LIBS': 'sundials_cvodes,sundials_nvecserial' + (
        ',sundials_sunlinsollapackdense,sundials_sunlinsollapackband' if _sun3 else ''
    ),
}
