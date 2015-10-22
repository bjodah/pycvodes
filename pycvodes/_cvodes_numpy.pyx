# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
cimport numpy as cnp
import numpy as np

from cvodes_numpy cimport PyCvodes

cnp.import_array()  # Numpy C-API initialization

cdef class Cvodes:

    cdef PyCvodes *thisptr

    def __cinit__(self, object f, object j, size_t ny, int ml=-1, int mu=-1):
        self.thisptr = new PyCvodes(<PyObject *>f, <PyObject *>j, ny, ml, mu)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                 double t0, double tend,
                 double atol, double rtol,
                 double hstart=0.0, int step_type_idx=1):
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        return self.thisptr.adaptive(<PyObject*>y0, t0, tend, atol,
                                     rtol, hstart, step_type_idx)

    def predefined(self, cnp.ndarray[cnp.float64_t, ndim=1] y0,
                   cnp.ndarray[cnp.float64_t, ndim=1] xout,
                   double dx0, double atol, double rtol,
                   int step_type_idx=8, double dx_max=0, double dx_min=0):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((xout.size, y0.size),
                                                                dtype=np.float64)
        if y0.size < self.thisptr.ny:
            raise ValueError("y0 too short")
        yout[0, :] = y0
        self.thisptr.predefined(<PyObject*>y0, <PyObject*>xout, <PyObject*>yout,
                                dx0, atol, rtol, step_type_idx, dx_max, dx_min)
        return yout

    def get_xout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] xout = np.empty(nsteps, dtype=np.float64)
        cdef size_t i
        for i in range(nsteps):
            xout[i] = self.thisptr.xout[i]
        return xout

    def get_yout(self, size_t nsteps):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] yout = np.empty((nsteps, self.thisptr.ny),
                                                                dtype=np.float64)
        cdef size_t i
        cdef size_t ny = self.thisptr.ny
        for i in range(nsteps):
            for j in range(ny):
                yout[i, j] = self.thisptr.yout[i*ny + j]
        return yout

    def get_info(self):
        return {'nrhs': self.thisptr.nrhs, 'njac': self.thisptr.njac}


steppers = ['adams', 'bdf']
requires_jac = ('bdf',)

def adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, method='bdf',
             lband=None, uband=None):
    cdef size_t nsteps
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = Cvodes(rhs, jac, len(y0),
                    -1 if lband is None else lband,
                    -1 if uband is None else uband)
    nsteps = integr.adaptive(np.array(y0, dtype=np.float64),
                             x0, xend, dx0, atol, rtol,
                             steppers.index(method))
    return integr.get_xout(nsteps), integr.get_yout(nsteps), integr.get_info()


def predefined(rhs, jac, y0, xout, dx0, atol, rtol, method='bdf',
               lband=None, uband=None):
    if method in requires_jac and jac is None:
        raise ValueError("Method requires explicit jacobian callback")
    integr = Cvodes(rhs, jac, len(y0),
                    -1 if lband is None else lband,
                    -1 if uband is None else uband)
    yout = integr.predefined(np.array(y0, dtype=np.float64),
                             np.array(xout, dtype=np.float64),
                             dx0, atol, rtol, steppers.index(method))
    return yout, integr.get_info()
