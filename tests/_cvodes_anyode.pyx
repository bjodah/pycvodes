# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# distutils: extra_compile_args = ['-std=c++14']

from libc.stdlib cimport malloc
from libcpp.vector cimport vector
from cvodes_cxx cimport Adams
from cvodes_anyode cimport simple_adaptive
cimport numpy as cnp
import numpy as np

cdef extern from "cvodes_cxx.hpp":
    ctypedef double realtype

cnp.import_array()  # Numpy C-API initialization

cdef extern from "testing_utils.hpp":
    cppclass Decay:
        Decay(realtype)

cdef class PyDecay:
    cdef Decay *thisptr

    def __cinit__(self, realtype k):
        self.thisptr = new Decay(k)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, realtype y0, realtype t):
        cdef:
            vector[int] root_indices
            realtype [:,::1] xyout_view
            cnp.npy_intp xyout_dims[2]
            int nout
            int td = 1
            realtype * xyout = <realtype *>malloc(2*sizeof(realtype))
        xyout[0] = 0
        xyout[1] = y0
        nout = simple_adaptive[Decay](&xyout, &td, self.thisptr, [1e-10], 1e-10, Adams, t, root_indices)
        xyout_dims[0] = nout + 1
        xyout_dims[1] = 2
        xyout_view = <realtype [:xyout_dims[0], :xyout_dims[1]]> xyout
        xyout_arr = np.asarray(xyout_view)
        xout = xyout_arr[:, 0]
        yout = xyout_arr[:, 1:]
        return xout, yout
