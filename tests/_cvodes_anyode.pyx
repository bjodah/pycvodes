# -*- coding: utf-8; mode: cython -*-
# distutils: language = c++
# distutils: extra_compile_args = ['-std=c++11']

from libcpp.vector cimport vector
from cvodes_cxx cimport Adams
from cvodes_anyode cimport simple_adaptive


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
        return simple_adaptive[Decay[double]](self.thisptr, [1e-10], 1e-10, Adams, &y0, 0.0, t, root_indices)
