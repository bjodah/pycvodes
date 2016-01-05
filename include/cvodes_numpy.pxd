# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "cvodes_numpy.hpp" namespace "cvodes_numpy":
    cdef cppclass PyCvodes:
        const size_t ny
        size_t nrhs, njac
        double time_cpu
        const int mlower, mupper, nroots
        vector[double] xout
        vector[double] yout
        vector[int] root_indices

        PyCvodes(PyObject*, PyObject*, PyObject*, size_t, int, int, int)
        size_t adaptive(PyObject*, double, double, double, double, int,
                        double, double, double, long int, int, int, int, bool) except +
        void predefined(PyObject*, PyObject*, PyObject*, double, double, int,
                        double, double, double, long int, int, int) except +
