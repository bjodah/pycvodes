# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp.vector cimport vector

cdef extern from "cvodes_numpy.hpp" namespace "cvodes_numpy":
    cdef cppclass PyCvodesOdeiv2:
        size_t ny
        vector[double] xout
        vector[double] yout

        PyCvodesOdeiv2(PyObject*, PyObject*, size_t)
        size_t adaptive(PyObject*, double, double, double, double, double, int) except +
        void predefined(PyObject*, PyObject*, PyObject*, double, double, double, int,
                        double, double) except +
