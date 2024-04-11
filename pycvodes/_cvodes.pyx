# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# cython: language_level=3str

import warnings
from cpython.object cimport PyObject
from libc.stdlib cimport malloc
from libcpp cimport bool
from libcpp.vector cimport vector
from collections.abc import Iterable

from cvodes_anyode_numpy cimport CvodesPyOdeSys
from cvodes_cxx cimport lmm_from_name, iter_type_from_name, \
     linear_solver_from_name, fpes as _fpes
from cvodes_anyode cimport simple_adaptive, simple_predefined

import numpy as np
cimport numpy as cnp

# Need these here rather than as imports so that
# typedefs are available at compile- (not run-) time.
# NOTE: base types of "int" and "float" are just
# appropriately-close standins as per Cython rules; will
# be replaced with the exact extern typedef at compile-time.
cdef extern from "cvodes_cxx.hpp":
     ctypedef double realtype
     ctypedef int indextype
     int PYCVODES_NO_KLU, PYCVODES_NO_LAPACK

cdef extern from "sundials_cxx.hpp" namespace "sundials_cxx":
    int version_major, version_minor, version_patch

cnp.import_array()  # Numpy C-API initialization

steppers = ('adams', 'bdf')
requires_jac = ('bdf',)
iterative_linsols = ('gmres', 'gmres_classic', 'bicgstab', 'tfqmr')
sundials_version = (version_major, version_minor, version_patch)

env = {
    "KLU": PYCVODES_NO_KLU != 1,
    "LAPACK": PYCVODES_NO_LAPACK != 1
}

fpes = {str(k.decode('utf-8')): v for k, v in dict(_fpes).items()}

# These need to be available as type objects at run type, in addition to the corresponding
# type tags (e.g. np.float64_t), which only exist at compile time and cannot be used with
# np.asarray(..., dtype=)
if sizeof(realtype) == sizeof(cnp.npy_double):
    dtype = np.float64
    env["REAL_TYPE"] = "double"
    env["SUNDIALS_PRECISION"] = "double"
elif sizeof(realtype) == sizeof(cnp.npy_float):
    dtype = np.float32
    env["REAL_TYPE"] = "float"
    env["SUNDIALS_PRECISION"] = "single"
elif sizeof(realtype) == sizeof(cnp.npy_longdouble):
    dtype = np.longdouble
    env["REAL_TYPE"] = "long double"
    env["SUNDIALS_PRECISION"] = "extended"
else:
    dtype = np.float64
    env["REAL_TYPE"] = "realtype"   # unclear

if sizeof(indextype) == sizeof(cnp.npy_int):
    env["INDEX_TYPE"] = "int"
else:
    env["INDEX_TYPE"] = "long int"


# signature in python methods should be able to accept any floating type regardless
# of what realtype is under the hood. scalars of type "floating" passed to the cython wrapper
# should be auto-cast to realtype when passed to C functions; any vectors/arrays
# will be manually cast below
ctypedef fused floating:
    cnp.float32_t
    cnp.float64_t
    cnp.longdouble_t

# cdef dict get_last_info(CvodesPyOdeSys * odesys, success=True):
#     info = {str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_int).items()}
#     info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_dbl).items()})
#     info.update({str(k.decode('utf-8')): np.array(v, dtype=np.float64) for k, v in dict(odesys.current_info.nfo_vecdbl).items()})
#     info.update({str(k.decode('utf-8')): np.array(v, dtype=int) for k, v in dict(odesys.current_info.nfo_vecint).items()})
#     info['nfev'] = odesys.nfev
#     info['njev'] = odesys.njev
#     info['success'] = success
#     return info

cdef _reshape_roots(cnp.ndarray roots, indextype ny):
    cdef cnp.ndarray out = roots.reshape((roots.size // (ny + 1), ny + 1))
    return out[:, 0], out[:, 1:]

def adaptive(rhs, jac, floating [:] yq0, floating x0, floating xend, atol,
             rtol, str method='bdf'
             , int nsteps=500, floating dx0=0.0, floating dx_min=0.0,
             floating dx_max=0.0, const unsigned nderiv=0, roots=None, int nroots=0, bool return_on_root=False
             ):
    cdef:
        int lband=-1, uband=-1, nquads=0
        indextype nyq = yq0.shape[yq0.ndim - 1]
        indextype ny = nyq - nquads
        bool with_jacobian = jac is not None
        bint with_jtimes = False
        CvodesPyOdeSys * odesys
        str iter_type="undecided"
        str linear_solver="default"
        int maxl=0,
        floating eps_lin=0.0
        int autorestart=0, return_on_error=False, record_rhs_xvals=False
        bool record_jac_xvals=False, record_order=False, record_fpe=False
        bool record_steps=False, autonomous_exprs=False
        int td = 500
        int tidx = 0
        realtype * xyqout = <realtype * const>malloc(td*(1 + ny*(nderiv+1) + nquads)*sizeof(realtype))
        realtype * ew_ele_out = NULL
        realtype [:, ::1] xyqout_view
        cnp.npy_intp xyqout_dims[2]
        int nout = 2
        int i

    dx_max_cb=None
    dx0cb=None

    xyqout[0] = <realtype> x0
    for i in range(ny):
        xyqout[1+i] = <realtype> yq0[i]

    for i in range(nquads):
        xyqout[1+ny*(nderiv+1)+i] = 0.0;
    nnz = -1
    odesys = new CvodesPyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *> None, <PyObject *>None,
                          <PyObject *>None, <PyObject *>None, lband, uband, nquads, nroots,
                          <PyObject *>dx0cb, <PyObject *>dx_max_cb, nnz)
    odesys.autonomous_exprs = autonomous_exprs
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    odesys.record_steps = record_steps

    try:
        xyqout_dims[0] = nout + 1
        xyqout_dims[1] = ny*(nderiv+1) + 1 + nquads
        xyqout_view = <realtype [:xyqout_dims[0], :xyqout_dims[1]]> xyqout
        xyqout_arr = np.asarray(xyqout_view)

        xout = xyqout_arr[:, 0]
        yout = xyqout_arr[:, 1:1+ny*(nderiv+1)]
        if return_on_error:
            success = xout[-1] == xend
        else:
            success = True

        info = {} #get_last_info(odesys, success)
        info['atol'] = [atol]
        info['rtol'] = rtol
        # if nroots > 0:
        #     info['root_indices'] = root_indices
        if nquads > 0:
            info['quads'] = xyqout_arr[:, 1+(1+nderiv)*ny:]
        yout_shape = (xout.size, ny) if nderiv == 0 else (xout.size, nderiv+1, ny)

        return xout, yout.reshape(yout_shape), info
    finally:
        del odesys

def predefined(rhs, jac,
               floating [::1] yq0,
               floating [::1] xout, atol,
               floating rtol, str method='bdf', int nsteps=500, floating dx0=0.0, floating dx_min=0.0,
               floating dx_max=0.0, quads=None, roots=None, cb_kwargs=None, int lband=-1, int uband=-1,
               int nquads=0, int nroots=0, str iter_type="undecided", str linear_solver="default", const int maxl=0,
               const floating eps_lin=0.0, const unsigned nderiv=0, bool return_on_root=False,
               int autorestart=0, bool return_on_error=False, bool record_rhs_xvals=False,
               bool record_jac_xvals=False, bool record_order=False, bool record_fpe=False,
               bool record_steps=False, dx0cb=None, dx_max_cb=None, bool autonomous_exprs=False,
               jtimes=None, bool ew_ele=False, indextype nnz=-1, const vector[realtype] constraints=[],
               long int max_num_steps_between_jac=0, bool stab_lim_det=False):
    return 0
