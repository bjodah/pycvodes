# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++

import warnings
from cpython.object cimport PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport numpy as cnp

from anyode_numpy cimport PyOdeSys
from cvodes_cxx cimport lmm_from_name, iter_type_from_name
from cvodes_anyode cimport simple_adaptive, simple_predefined

import numpy as np

cnp.import_array()  # Numpy C-API initialization

steppers = ('adams', 'bdf')
requires_jac = ('bdf',)

iter_types = {'default': 0, 'functional': 1, 'newton': 2}  # grep "define CV_FUNCTIONAL" cvodes.h
linear_solvers = {'default': 0, 'dense': 1, 'banded': 2, 'gmres': 10, 'gmres_classic': 11, 'bicgstab': 20, 'tfqmr': 30}


cdef dict get_last_info(PyOdeSys * odesys, success=True):
    info = {str(k.decode('utf-8')): v for k, v in dict(odesys.last_integration_info).items()}
    info.update({str(k.decode('utf-8')): v for k, v in dict(odesys.last_integration_info_dbl).items()})
    info['nfev'] = odesys.nfev
    info['njev'] = odesys.njev
    info['success'] = success
    return info

cdef _reshape_roots(cnp.ndarray[cnp.float64_t, ndim=1] roots, int ny):
    cdef cnp.ndarray[cnp.float64_t, ndim=2] out = roots.reshape((roots.size / (ny + 1), ny + 1))
    return out[:, 0], out[:, 1:]


def adaptive(rhs, jac, cnp.ndarray[cnp.float64_t, mode='c'] y0, double x0, double xend, double atol,
             double rtol, str method='bdf', int nsteps=500, double dx0=0.0, double dx_min=0.0,
             double dx_max=0.0, roots=None, cb_kwargs=None, int lband=-1, int uband=-1, int nroots=0,
             str iter_type="undecided", int linear_solver=0, const int maxl=0,
             const double eps_lin=0.0, const unsigned nderiv=0, bool return_on_root=False,
             int autorestart=0, bool return_on_error=False):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        bool with_jacobian = jac is not None
        PyOdeSys * odesys
        vector[int] root_indices

    if method.lower() in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    if np.isinf([x0, xend]).any(): raise ValueError("+/-Inf found in x0/xend")
    if np.isnan([x0, xend]).any(): raise ValueError("NaN found in x0/xend")
    if np.isinf(y0).any(): raise ValueError("+/-Inf found in y0")
    if np.isnan(y0).any(): raise ValueError("NaN found in y0")

    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *>roots,
                          <PyObject *>cb_kwargs, lband, uband, nroots)
    try:
        xout, yout = map(np.asarray, simple_adaptive[PyOdeSys](
            odesys, [atol], rtol, lmm_from_name(method.lower().encode('UTF-8')),
            &y0[0], x0, xend, root_indices, nsteps, dx0, dx_min, dx_max, with_jacobian,
            iter_type_from_name(iter_type.lower().encode('UTF-8')),
            linear_solver, maxl, eps_lin, nderiv, return_on_root, autorestart, return_on_error))
        info = get_last_info(odesys, False if return_on_error and xout[-1] != xend else True)
        if nroots > 0:
            info['root_indices'] = root_indices
        yout_shape = (xout.size, ny) if nderiv == 0 else (xout.size, nderiv+1, ny)
        return xout, yout.reshape(yout_shape), info
    finally:
        del odesys


def predefined(rhs, jac,
               cnp.ndarray[cnp.float64_t, mode='c'] y0,
               cnp.ndarray[cnp.float64_t, ndim=1] xout,
               double atol, double rtol, str method='bdf', int nsteps=500, double dx0=0.0,
               double dx_min=0.0, double dx_max=0.0, roots=None, cb_kwargs=None, int lband=-1, int uband=-1, int nroots=0,
               str iter_type="undecided", int linear_solver=0, const int maxl=0, const double eps_lin=0.0,
               const unsigned nderiv=0, bool return_on_root=False, int autorestart=0, bool return_on_error=False):
    cdef:
        int ny = y0.shape[y0.ndim - 1]
        cnp.ndarray[cnp.float64_t, ndim=3] yout = np.empty((xout.size, nderiv+1, ny))
        bool with_jacobian = jac is not None
        int nreached
        PyOdeSys * odesys
        vector[int] root_indices
        vector[double] roots_output

    if method.lower() in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    if np.isinf(xout).any(): raise ValueError("+/-Inf found in xout")
    if np.isnan(xout).any(): raise ValueError("NaN found in xout")
    if np.isinf(y0).any(): raise ValueError("+/-Inf found in y0")
    if np.isnan(y0).any(): raise ValueError("NaN found in y0")
    odesys = new PyOdeSys(ny, <PyObject *>rhs, <PyObject *>jac, <PyObject *>roots,
                          <PyObject *>cb_kwargs, lband, uband, nroots)
    try:
        nreached = simple_predefined[PyOdeSys](
            odesys, [atol], rtol, lmm_from_name(method.lower().encode('UTF-8')), &y0[0],
            xout.size, &xout[0], <double *>yout.data, root_indices, roots_output, nsteps,
            dx0, dx_min, dx_max, with_jacobian, iter_type_from_name(iter_type.lower().encode('UTF-8')),
            linear_solver, maxl, eps_lin, nderiv, autorestart, return_on_error)
        info = get_last_info(odesys, success=False if return_on_error and nreached < xout.size else True)
        info['nreached'] = nreached
        if nroots > 0:
            info['root_indices'] = root_indices
            info['roots_output'] = _reshape_roots(np.asarray(roots_output), ny)
        return yout.reshape((xout.size, ny)) if nderiv == 0 else yout, info
    finally:
        del odesys
