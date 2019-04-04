# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# cython: language_level=3str

import warnings
from cpython.object cimport PyObject
from libc.stdlib cimport malloc
from libcpp cimport bool
from libcpp.vector cimport vector
from collections import Iterable

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

cdef extern from "sundials_cxx.hpp" namespace "sundials_cxx":
    int version_major, version_minor, version_patch

cnp.import_array()  # Numpy C-API initialization

steppers = ('adams', 'bdf')
requires_jac = ('bdf',)
iterative_linsols = ('gmres', 'gmres_classic', 'bicgstab', 'tfqmr')
sundials_version = (version_major, version_minor, version_patch)


fpes = {str(k.decode('utf-8')): v for k, v in dict(_fpes).items()}

# These need to be available as type objects at run type, in addition to the corresponding
# type tags (e.g. np.float64_t), which only exist at compile time and cannot be used with
# np.asarray(..., dtype=)
if sizeof(realtype) == sizeof(cnp.npy_double):
    dtype = np.float64
elif sizeof(realtype) == sizeof(cnp.npy_float):
    dtype = np.float32
elif sizeof(realtype) == sizeof(cnp.npy_longdouble):
    dtype = np.longdouble
else:
    dtype = np.float64

# signature in python methods should be able to accept any floating type regardless
# of what realtype is under the hood. scalars of type "floating" passed to the cython wrapper
# should be auto-cast to realtype when passed to C functions; any vectors/arrays
# will be manually cast below
ctypedef fused floating:
    cnp.float32_t
    cnp.float64_t
    cnp.longdouble_t

cdef dict get_last_info(CvodesPyOdeSys * odesys, success=True):
    info = {str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_int).items()}
    info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_dbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=np.float64) for k, v in dict(odesys.current_info.nfo_vecdbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=int) for k, v in dict(odesys.current_info.nfo_vecint).items()})
    info['nfev'] = odesys.nfev
    info['njev'] = odesys.njev
    info['success'] = success
    return info

cdef _reshape_roots(cnp.ndarray roots, indextype ny):
    cdef cnp.ndarray out = roots.reshape((roots.size // (ny + 1), ny + 1))
    return out[:, 0], out[:, 1:]

def adaptive(rhs, jac, floating [:] yq0, floating x0, floating xend, atol,
             rtol, str method='bdf', int nsteps=500, floating dx0=0.0, floating dx_min=0.0,
             floating dx_max=0.0, quads=None, roots=None, cb_kwargs=None, int lband=-1, int uband=-1,
             int nquads=0, int nroots=0,
             str iter_type="undecided", str linear_solver="default", const int maxl=0,
             const floating eps_lin=0.0, const unsigned nderiv=0, bool return_on_root=False,
             int autorestart=0, bool return_on_error=False, bool record_rhs_xvals=False,
             bool record_jac_xvals=False, bool record_order=False, bool record_fpe=False,
             bool record_steps=False, dx0cb=None, dx_max_cb=None, bool autonomous_exprs=False,
             int nprealloc=500, jtimes=None, bool ew_ele=False, indextype nnz=-1, vector[realtype] constraints=[]):
    cdef:
        indextype nyq = yq0.shape[yq0.ndim - 1]
        indextype ny = nyq - nquads
        bool with_jacobian = jac is not None
        bool with_jtimes = jtimes is not None
        CvodesPyOdeSys * odesys
        vector[int] root_indices
        vector[realtype] atol_vec
        int td = nprealloc
        int tidx = 0
        realtype * xyqout = <realtype * const>malloc(td*(1 + ny*(nderiv+1) + nquads)*sizeof(realtype))
        realtype * ew_ele_out = NULL
        realtype [:, ::1] xyqout_view
        realtype [:, :, ::1] ew_ele_view
        cnp.npy_intp xyqout_dims[2]
        cnp.npy_intp ew_ele_dims[3]
        int nout
        int i

    if nprealloc < 1:
        raise ValueError("Need nprealloc > 0")
    if isinstance(atol, Iterable):
        for at in atol:
            atol_vec.push_back(<realtype> at)
    else:
        atol_vec.push_back(<realtype> atol)
    rhs(0, yq0[..., :ny], np.empty(ny))  # fail early if rhs does not work

    if method.lower() in requires_jac and not with_jacobian:
        warnings.warn("No full jacobian provided; disabling default preconditioning.")
        if linear_solver.lower() not in iterative_linsols:
            warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
        elif not with_jtimes:
            warnings.warn("Method requires jacobian or jacobian-vector product, no callback provided: using finite differences (may be inaccurate).")

    if np.isinf([x0, xend]).any(): raise ValueError("+/-Inf found in x0/xend")
    if np.isnan([x0, xend]).any(): raise ValueError("NaN found in x0/xend")
    if np.isinf(yq0).any(): raise ValueError("+/-Inf found in yq0")
    if np.isnan(yq0).any(): raise ValueError("NaN found in yq0")
    xyqout[0] = <realtype> x0
    for i in range(ny):
        xyqout[1+i] = <realtype> yq0[i]
    if ew_ele:
        ew_ele_out = <realtype*>malloc(2*td*ny*sizeof(realtype))
        for i in range(ny):
            ew_ele_out[i] = 0.0

    for i in range(nquads):
        xyqout[1+ny*(nderiv+1)+i] = 0.0;

    odesys = new CvodesPyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *> jtimes, <PyObject *>quads,
                          <PyObject *>roots, <PyObject *>cb_kwargs, lband, uband, nquads, nroots,
                          <PyObject *>dx0cb, <PyObject *>dx_max_cb, nnz)
    odesys.autonomous_exprs = autonomous_exprs
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    odesys.record_steps = record_steps

    try:
        nout = simple_adaptive[CvodesPyOdeSys](
            &xyqout, &td, odesys, atol_vec, rtol, lmm_from_name(method.lower().encode('UTF-8')),
            xend, root_indices, nsteps, dx0, dx_min, dx_max, with_jacobian,
            iter_type_from_name(iter_type.lower().encode('UTF-8')),
            linear_solver_from_name(linear_solver.lower().encode('UTF-8')),
            maxl, eps_lin, nderiv, return_on_root, autorestart, return_on_error, with_jtimes,
            tidx, &ew_ele_out if ew_ele else NULL, constraints)

        xyqout_dims[0] = nout + 1
        xyqout_dims[1] = ny*(nderiv+1) + 1 + nquads
        xyqout_view = <realtype [:xyqout_dims[0], :xyqout_dims[1]]> xyqout
        xyqout_arr = np.asarray(xyqout_view)

        xout = xyqout_arr[:, 0]
        yout = xyqout_arr[:, 1:1+ny*(nderiv+1)]
        if return_on_error:
            if return_on_root and root_indices.size() > 0:
                success = root_indices[root_indices.size() - 1] == len(xout) - 1
            else:
                success = xout[-1] == xend
        else:
            success = True

        info = get_last_info(odesys, success)
        info['atol'] = atol_vec
        info['rtol'] = rtol
        if nroots > 0:
            info['root_indices'] = root_indices
        if nquads > 0:
            info['quads'] = xyqout_arr[:, 1+(1+nderiv)*ny:]
        yout_shape = (xout.size, ny) if nderiv == 0 else (xout.size, nderiv+1, ny)

        if ew_ele:
            ew_ele_dims[0] = nout + 1
            ew_ele_dims[1] = 2
            ew_ele_dims[2] = ny
            ew_ele_view = <realtype [:ew_ele_dims[0], :ew_ele_dims[1], :ew_ele_dims[2]]> ew_ele_out
            ew_ele_arr = np.asarray(ew_ele_view)
            info['ew_ele'] = ew_ele_arr
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
               jtimes=None, bool ew_ele=False, indextype nnz=-1, const vector[realtype] constraints=[]):
    cdef:
        indextype nyq = yq0.shape[yq0.ndim - 1]
        indextype ny = nyq - nquads
        realtype * yqout = <realtype *>malloc(xout.size*(nderiv + 1)*nyq*sizeof(realtype))
        realtype * ew_ele_out = NULL
        bool with_jacobian = jac is not None
        bool with_jtimes = jtimes is not None
        int nreached
        CvodesPyOdeSys * odesys
        vector[int] root_indices
        vector[realtype] roots_output
        vector[realtype] atol_vec
        realtype [:, :, ::1] yqout_view
        realtype [:, :, ::1] ew_ele_view
        cnp.npy_intp yqout_dims[3]
        cnp.npy_intp ew_ele_dims[3]

    rhs(0, yq0[..., :ny], np.empty(ny))  # fail early if rhs does not work

    if isinstance(atol, Iterable):
        for at in atol:
            atol_vec.push_back(<realtype> at)
    else:
        atol_vec.push_back(<realtype> atol)

    if ew_ele:
        ew_ele_out = <realtype*>malloc(xout.size*2*ny*sizeof(realtype))

    if method.lower() in requires_jac and not with_jacobian:
        warnings.warn("No full jacobian provided; disabling default preconditioning.")
        if linear_solver.lower() not in iterative_linsols:
            warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
        elif not with_jtimes:
            warnings.warn("Method requires jacobian or jacobian-vector product, no callback provided: using finite differences (may be inaccurate).")

    if np.isinf(xout).any(): raise ValueError("+/-Inf found in xout")
    if np.isnan(xout).any(): raise ValueError("NaN found in xout")
    if np.isinf(yq0).any(): raise ValueError("+/-Inf found in yq0")
    if np.isnan(yq0).any(): raise ValueError("NaN found in yq0")

    odesys = new CvodesPyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *> jtimes, <PyObject *>quads,
                          <PyObject *>roots, <PyObject *>cb_kwargs, lband, uband, nquads, nroots,
                          <PyObject *>dx0cb, <PyObject *>dx_max_cb, nnz)
    odesys.autonomous_exprs = autonomous_exprs
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    odesys.record_steps = record_steps

    # this process is necessary since maybe the input type floating != realtype
    cdef cnp.ndarray[realtype, ndim=1, mode='c'] yq0_arr = np.asarray(yq0, dtype=dtype)
    cdef cnp.ndarray[realtype, ndim=1, mode='c'] xout_arr = np.asarray(xout, dtype=dtype)

    try:
        nreached = simple_predefined[CvodesPyOdeSys](
            odesys, atol_vec, rtol, lmm_from_name(method.lower().encode('UTF-8')), &yq0_arr[0],
            xout.size, &xout_arr[0], yqout, root_indices, roots_output, nsteps,
            dx0, dx_min, dx_max, with_jacobian, iter_type_from_name(iter_type.lower().encode('UTF-8')),
            linear_solver_from_name(linear_solver.lower().encode('UTF-8')),
            maxl, eps_lin, nderiv, autorestart, return_on_error, with_jtimes,
            ew_ele_out if ew_ele else NULL, constraints)

        yqout_dims[0] = xout.size
        yqout_dims[1] = nderiv + 1
        yqout_dims[2] = nyq
        yqout_view = <realtype [:yqout_dims[0], :yqout_dims[1], :yqout_dims[2]:1]> yqout
        yqout_arr = np.asarray(yqout_view)

        info = get_last_info(odesys, success=False if return_on_error and nreached < xout.size else True)
        info['nreached'] = nreached
        info['atol'] = atol_vec
        info['rtol'] = rtol
        if nquads > 0:
            info['quads'] = yqout_arr[:, 0, ny:]
        if nroots > 0:
            info['root_indices'] = root_indices
            info['roots_output'] = _reshape_roots(np.asarray(roots_output), ny)
        if ew_ele:
            ew_ele_dims[0] = xout.size
            ew_ele_dims[1] = 2
            ew_ele_dims[2] = ny
            ew_ele_view = <realtype [:ew_ele_dims[0], :ew_ele_dims[1], :ew_ele_dims[2]]> ew_ele_out
            ew_ele_arr = np.asarray(ew_ele_view)
            info['ew_ele'] = ew_ele_arr
        yout = yqout_arr[:, :, :ny]
        return yout.reshape((xout.size, ny)) if nderiv == 0 else yout, info
    finally:
        del odesys
