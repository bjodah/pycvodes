# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# distutils: extra_compile_args = ['-std=c++11']

from libc.stdlib cimport malloc
from libcpp.vector cimport vector
from cvodes_cxx cimport Adams
from cvodes_anyode cimport simple_adaptive
cimport numpy as cnp

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(cnp.ndarray arr, int flags)

cnp.import_array()  # Numpy C-API initialization

cdef extern from "testing_utils.hpp":
    cppclass Decay[T]:
        Decay(T)


cdef class PyDecay:
    cdef Decay[double] *thisptr

    def __cinit__(self, double k):
        self.thisptr = new Decay[double](k)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, double y0, double t):
        cdef:
            vector[int] root_indices
            cnp.ndarray[cnp.float64_t, ndim=2] xyout_arr
            cnp.npy_intp xyout_dims[2]
            int nout
            int td = 1
            double * xyout = <double *>malloc(2*sizeof(double))
        xyout[0] = 0
        xyout[1] = y0
        nout = simple_adaptive[Decay[double]](&xyout, &td, self.thisptr, [1e-10], 1e-10, Adams, t, root_indices)
        xyout_dims[0] = nout + 1
        xyout_dims[1] = 2
        xyout_arr = cnp.PyArray_SimpleNewFromData(
            2, xyout_dims, cnp.NPY_DOUBLE, <void *>xyout)
        PyArray_ENABLEFLAGS(xyout_arr, cnp.NPY_OWNDATA)
        xout = xyout_arr[:, 0]
        yout = xyout_arr[:, 1:]
        return xout, yout
