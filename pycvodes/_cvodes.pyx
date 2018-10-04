# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++

import warnings
from cpython.object cimport PyObject
from libc.stdlib cimport malloc
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport numpy as cnp

from anyode_numpy cimport PyOdeSys
from cvodes_cxx cimport lmm_from_name, iter_type_from_name, fpes as _fpes
from cvodes_anyode cimport simple_adaptive, simple_predefined

import numpy as np

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(cnp.ndarray arr, int flags)

cdef extern from "sundials_cxx.hpp" namespace "sundials_cxx":
    int version_major, version_minor, version_patch

cnp.import_array()  # Numpy C-API initialization

steppers = ('adams', 'bdf')
requires_jac = ('bdf',)
sundials_version = (version_major, version_minor, version_patch)

iter_types = {'default': 0, 'functional': 1, 'newton': 2}  # grep "define CV_FUNCTIONAL" cvodes.h
linear_solvers = {'default': 0, 'dense': 1, 'banded': 2, 'gmres': 10, 'gmres_classic': 11, 'bicgstab': 20, 'tfqmr': 30}

fpes = {str(k.decode('utf-8')): v for k, v in dict(_fpes).items()}

cdef dict get_last_info(PyOdeSys * odesys, success=True):
    info = {str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_int).items()}
    info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.current_info.nfo_dbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=np.float64) for k, v in dict(odesys.current_info.nfo_vecdbl).items()})
    info.update({str(k.decode('utf-8')): np.array(v, dtype=int) for k, v in dict(odesys.current_info.nfo_vecint).items()})
    info['nfev'] = odesys.nfev
    info['njev'] = odesys.njev
    info['success'] = success
    return info

cdef _reshape_roots(cnp.ndarray[cnp.float64_t, ndim=1] roots, int ny):
    cdef cnp.ndarray[cnp.float64_t, ndim=2] out = roots.reshape((roots.size // (ny + 1), ny + 1))
    return out[:, 0], out[:, 1:]


def adaptive(rhs, jac, cnp.ndarray[cnp.float64_t, mode='c'] yq0, double x0, double xend, atol,
             double rtol, str method='bdf', int nsteps=500, double dx0=0.0, double dx_min=0.0,
             double dx_max=0.0, quads=None, roots=None, cb_kwargs=None, int lband=-1, int uband=-1,
             int nquads=0, int nroots=0,
             str iter_type="undecided", int linear_solver=0, const int maxl=0,
             const double eps_lin=0.0, const unsigned nderiv=0, bool return_on_root=False,
             int autorestart=0, bool return_on_error=False, bool record_rhs_xvals=False,
             bool record_jac_xvals=False, bool record_order=False, bool record_fpe=False,
             bool record_steps=False, dx0cb=None, dx_max_cb=None, bool autonomous_exprs=False, int nprealloc=500,
             bool with_jtimes=False, bool ew_ele=False, vector[double] constraints=[]):
    cdef:
        int nyq = yq0.shape[yq0.ndim - 1]
        int ny = nyq - nquads
        bool with_jacobian = jac is not None
        PyOdeSys * odesys
        vector[int] root_indices
        vector[double] atol_vec
        int td = nprealloc
        int tidx = 0
        double * xyqout = <double*>malloc(td*(1 + ny*(nderiv+1) + nquads)*sizeof(double))
        double * ew_ele_out
        cnp.ndarray[cnp.float64_t, ndim=2] xyqout_arr
        cnp.ndarray[cnp.float64_t, ndim=3] ew_ele_arr
        cnp.npy_intp xyqout_dims[2]
        cnp.npy_intp ew_ele_dims[3]
        int nout
    if nprealloc < 1:
        raise ValueError("Need nprealloc > 0")
    if isinstance(atol, float):
        atol_vec.push_back(atol)
    else:
        for at in atol:
            atol_vec.push_back(at)
    rhs(0, yq0[..., :ny], np.empty(ny))  # fail early if rhs does not work


    if method.lower() in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    if np.isinf([x0, xend]).any(): raise ValueError("+/-Inf found in x0/xend")
    if np.isnan([x0, xend]).any(): raise ValueError("NaN found in x0/xend")
    if np.isinf(yq0).any(): raise ValueError("+/-Inf found in yq0")
    if np.isnan(yq0).any(): raise ValueError("NaN found in yq0")
    xyqout[0] = x0
    for i in range(ny):
        xyqout[1+i] = yq0[i]
    if ew_ele:
        ew_ele_out = <double*>malloc(2*td*ny*sizeof(double))
        for i in range(ny):
            ew_ele_out[i] = 0.0

    for i in range(nquads):
        xyqout[1+ny*(nderiv+1)+i] = 0.0;

    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *>quads, <PyObject *>roots,
                          <PyObject *>cb_kwargs, lband, uband, nquads, nroots, <PyObject *>dx0cb,
                          <PyObject *>dx_max_cb)
    odesys.autonomous_exprs = autonomous_exprs
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    odesys.record_steps = record_steps

    try:
        nout = simple_adaptive[PyOdeSys](
            &xyqout, &td, odesys, atol_vec, rtol, lmm_from_name(method.lower().encode('UTF-8')),
            xend, root_indices, nsteps, dx0, dx_min, dx_max, with_jacobian,
            iter_type_from_name(iter_type.lower().encode('UTF-8')), linear_solver, maxl,
            eps_lin, nderiv, return_on_root, autorestart, return_on_error, with_jtimes,
            tidx, &ew_ele_out if ew_ele else NULL, constraints)

        xyqout_dims[0] = nout + 1
        xyqout_dims[1] = ny*(nderiv+1) + 1 + nquads
        xyqout_arr = cnp.PyArray_SimpleNewFromData(2, xyqout_dims, cnp.NPY_DOUBLE, <void *>xyqout)
        PyArray_ENABLEFLAGS(xyqout_arr, cnp.NPY_OWNDATA)

        xout = xyqout_arr[:, 0]
        yout = xyqout_arr[:, 1:1+ny*(nderiv+1)]
        if return_on_error:
            if return_on_root and root_indices[root_indices.size() - 1] == len(xout) - 1:
                success = True
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
            ew_ele_arr = cnp.PyArray_SimpleNewFromData(3, ew_ele_dims, cnp.NPY_DOUBLE, <void *>ew_ele_out)
            PyArray_ENABLEFLAGS(ew_ele_arr, cnp.NPY_OWNDATA)
            info['ew_ele'] = ew_ele_arr
        return xout, yout.reshape(yout_shape), info
    finally:
        del odesys


def predefined(rhs, jac,
               cnp.ndarray[cnp.float64_t, mode='c'] yq0,
               cnp.ndarray[cnp.float64_t, ndim=1] xout, atol,
               double rtol, str method='bdf', int nsteps=500, double dx0=0.0, double dx_min=0.0,
               double dx_max=0.0, quads=None, roots=None, cb_kwargs=None, int lband=-1, int uband=-1,
               int nquads=0, int nroots=0, str iter_type="undecided", int linear_solver=0, const int maxl=0,
               const double eps_lin=0.0, const unsigned nderiv=0, bool return_on_root=False,
               int autorestart=0, bool return_on_error=False, bool record_rhs_xvals=False,
               bool record_jac_xvals=False, bool record_order=False, bool record_fpe=False,
               bool record_steps=False, dx0cb=None, dx_max_cb=None, bool autonomous_exprs=False,
               bool with_jtimes=False, bool ew_ele=False, vector[double] constraints=[]):
    cdef:
        int nyq = yq0.shape[yq0.ndim - 1]
        int ny = nyq - nquads
        cnp.ndarray[cnp.float64_t, ndim=3] yqout = np.empty((xout.size, nderiv+1, nyq))
        cnp.ndarray[cnp.float64_t, ndim=3] ew_ele_arr = np.empty((xout.size, 2, ny))
        bool with_jacobian = jac is not None
        int nreached
        PyOdeSys * odesys
        vector[int] root_indices
        vector[double] roots_output
        vector[double] atol_vec
    rhs(0, yq0[..., :ny], np.empty(ny))  # fail early if rhs does not work


    if isinstance(atol, float):
        atol_vec.push_back(atol)
    else:
        for at in atol:
            atol_vec.push_back(at)

    if method.lower() in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    if np.isinf(xout).any(): raise ValueError("+/-Inf found in xout")
    if np.isnan(xout).any(): raise ValueError("NaN found in xout")
    if np.isinf(yq0).any(): raise ValueError("+/-Inf found in yq0")
    if np.isnan(yq0).any(): raise ValueError("NaN found in yq0")
    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *>quads, <PyObject *>roots,
                          <PyObject *>cb_kwargs, lband, uband, nquads, nroots, <PyObject *>dx0cb, <PyObject *>dx_max_cb)
    odesys.autonomous_exprs = autonomous_exprs
    odesys.record_rhs_xvals = record_rhs_xvals
    odesys.record_jac_xvals = record_jac_xvals
    odesys.record_order = record_order
    odesys.record_fpe = record_fpe
    odesys.record_steps = record_steps

    try:
        nreached = simple_predefined[PyOdeSys](
            odesys, atol_vec, rtol, lmm_from_name(method.lower().encode('UTF-8')), &yq0[0],
            xout.size, &xout[0], <double *>yqout.data, root_indices, roots_output, nsteps,
            dx0, dx_min, dx_max, with_jacobian, iter_type_from_name(iter_type.lower().encode('UTF-8')),
            linear_solver, maxl, eps_lin, nderiv, autorestart, return_on_error, with_jtimes,
            <double *>ew_ele_arr.data if ew_ele else NULL, constraints)
        info = get_last_info(odesys, success=False if return_on_error and nreached < xout.size else True)
        info['nreached'] = nreached
        info['atol'] = atol_vec
        info['rtol'] = rtol
        if nquads > 0:
            info['quads'] = yqout[:, 0, ny:]
        if nroots > 0:
            info['root_indices'] = root_indices
            info['roots_output'] = _reshape_roots(np.asarray(roots_output), ny)
        if ew_ele:
            info['ew_ele'] = ew_ele_arr
        yout = yqout[:, :, :ny]
        return yout.reshape((xout.size, ny)) if nderiv == 0 else yout, info
    finally:
        del odesys
