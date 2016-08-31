# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++

import warnings
from cpython.object cimport PyObject
from libcpp cimport bool
cimport numpy as cnp
import numpy as np

from cvodes_numpy cimport PyCvodes

cnp.import_array()  # Numpy C-API initialization

steppers = {'adams': 1, 'bdf': 2}  # grep "define CV_ADAMS" cvodes.h
requires_jac = ('bdf',)
iter_types = {'default': 0, 'functional': 1, 'newton': 2}  # grep "define CV_FUNCTIONAL" cvodes.h
linear_solvers = {'default': 0, 'dense': 1, 'banded': 2, 'gmres': 10, 'gmres_classic': 11, 'bicgstab': 20, 'tfqmr': 30}


cdef class Cvodes:

    cdef PyCvodes *thisptr
    cdef bool success

    def __cinit__(self, object f, object j, size_t ny, object roots=None,
                  int ml=-1, int mu=-1, int nroots=0):
        self.thisptr = new PyCvodes(<PyObject *>f, <PyObject *>j, <PyObject *>roots,
                                    ny, ml, mu, nroots)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                 double t0, double tend,
                 double atol, double rtol,
                 int iter_type, int linear_solver, int step_type_idx=1,
                 double dx0=.0, double dx_min=.0, double dx_max=.0, long int mxsteps=0,
                 int nderiv=0, bool return_on_root=False):
        cdef:
            int maxl = 0
            double epslin = 0.0
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        self.success = True
        try:
            return self.thisptr.adaptive(
                <PyObject*>y0, t0, tend, atol, rtol, step_type_idx, dx0, dx_min, dx_max, mxsteps,
                iter_type, linear_solver, maxl, epslin, nderiv, return_on_root)
        except:
            self.success = False
            raise

    def predefined(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                   cnp.ndarray[cnp.float64_t, ndim=1] xout,
                   double atol, double rtol,
                   int iter_type, int linear_solver, int step_type_idx=8,
                   double dx0=.0, double dx_min=.0, double dx_max=.0, long int mxsteps=0,
                   int nderiv=0):
        cdef:
            int maxl = 0
            double epslin = 0.0
            cnp.ndarray[cnp.float64_t, ndim=3] yout = np.empty((xout.size, nderiv+1, y0.size),
                                                               dtype=np.float64)
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        yout[0, :] = y0
        self.success = True
        try:
            self.thisptr.predefined(
                <PyObject*>y0, <PyObject*>xout, <PyObject*>yout, atol, rtol, step_type_idx, dx0, dx_min, dx_max,
                mxsteps, iter_type, linear_solver, maxl, epslin, nderiv)
        except:
            self.success = False
            raise
        return yout.reshape((xout.size, y0.size)) if nderiv == 0 else yout

    def get_xout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] xout = np.empty(nsteps, dtype=np.float64)
        cdef size_t i
        for i in range(nsteps):
            xout[i] = self.thisptr.xout[i]
        return xout

    def get_yout(self, size_t nsteps, int nderiv=0):
        cdef cnp.ndarray[cnp.float64_t, ndim=3] yout = np.empty((nsteps, nderiv+1, self.thisptr.ny),
                                                                dtype=np.float64)
        cdef size_t i
        cdef size_t ny = self.thisptr.ny
        for i in range(nsteps):
            for j in range(nderiv+1):
                for k in range(ny):
                    yout[i, j, k] = self.thisptr.yout[i*ny*(nderiv+1) + j*ny + k]
        return yout.reshape((nsteps, self.thisptr.ny)) if nderiv == 0 else yout

    def get_root_indices(self):
        return self.thisptr.root_indices

    def get_roots_output(self):
        cdef:
            int i, j
            int ny_plus_1 = self.thisptr.ny + 1
        if self.thisptr.roots_output.size() % ny_plus_1 != 0:
            raise ValueError("roots_output of incompatible dimension.")
        shape = (self.thisptr.roots_output.size() / ny_plus_1, ny_plus_1)
        cdef cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] arr = np.empty(shape)
        for i in range(self.thisptr.roots_output.size() / ny_plus_1):
            for j in range(ny_plus_1):
                arr[i, j] = self.thisptr.roots_output[i*ny_plus_1 + j]
        return arr[:, 0], arr[:, 1:]

    def get_info(self):
        info = {'nfev': self.thisptr.nfev,
                'njev': self.thisptr.njev,
                'time_cpu': self.thisptr.time_cpu,
                'time_wall': self.thisptr.time_wall,
                'success': self.success}
        if self.thisptr.nroots > 0:
            info['root_indices'] = self.get_root_indices()
            info['roots_output'] = self.get_roots_output()
        info.update(self.last_integration_info)
        return info

    property last_integration_info:
        def __get__(self):
            return self.thisptr.last_integration_info


def adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, dx_min=0.0, dx_max=0.0, mxsteps=0, nderiv=0, method='bdf',
             iter_type='default', linear_solver='default', lband=None, uband=None, roots=None, nroots=0,
             return_on_root=False):
    cdef size_t nsteps
    if method in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    integr = Cvodes(rhs, jac, len(y0), roots,
                    -1 if lband is None else lband,
                    -1 if uband is None else uband,
                    nroots)
    nsteps = integr.adaptive(np.array(y0, dtype=np.float64), x0, xend, atol, rtol,
                             iter_types.get(iter_type, 0), linear_solvers[linear_solver], steppers[method],
                             dx0, dx_min, dx_max, mxsteps, nderiv, return_on_root)
    return integr.get_xout(nsteps), integr.get_yout(nsteps, nderiv), integr.get_info()


def predefined(rhs, jac, y0, xout, dx0, atol, rtol, dx_min=0.0, dx_max=0.0, mxsteps=0, nderiv=0, method='bdf',
               iter_type='default', linear_solver='default', lband=None, uband=None, roots=None, nroots=0):
    if method in requires_jac and jac is None:
        warnings.warn("Method requires jacobian, no callback provided: using finite differences (may be inaccurate).")
    integr = Cvodes(rhs, jac, len(y0), roots,
                    -1 if lband is None else lband,
                    -1 if uband is None else uband,
                    nroots)
    yout = integr.predefined(np.array(y0, dtype=np.float64), np.array(xout, dtype=np.float64), atol, rtol,
                             iter_types.get(iter_type, 0), linear_solvers[linear_solver], steppers[method],
                             dx0, dx_min, dx_max, mxsteps, nderiv)
    return yout, integr.get_info()
