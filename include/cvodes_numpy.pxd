# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector

cdef extern from "cvodes_numpy.hpp" namespace "cvodes_numpy":
    cdef cppclass PyCvodes:
        const size_t ny
        size_t nrhs, njac
        int mlower, mupper
        vector[double] xout
        vector[double] yout

        PyCvodes(PyObject*, PyObject*, size_t, int, int)
        size_t adaptive(PyObject*, double, double, double, double, double, int) except +
        void predefined(PyObject*, PyObject*, PyObject*, double, double, double, int,
                        double, double) except +
