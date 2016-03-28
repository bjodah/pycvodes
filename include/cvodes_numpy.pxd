# -*- coding: utf-8; mode: cython -*-

from cpython.object cimport PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map

cdef extern from "cvodes_numpy.hpp" namespace "cvodes_numpy":
    cdef cppclass PyCvodes:
        const size_t ny
        size_t nfev, njev
        double time_cpu, time_wall
        const int mlower, mupper, nroots
        vector[double] xout
        vector[double] yout
        vector[int] root_indices
        vector[int] roots_output
        unordered_map[string, int] last_integration_info

        PyCvodes(PyObject*, PyObject*, PyObject*, size_t, int, int, int)
        size_t adaptive(PyObject*, double, double, double, double, int,
                        double, double, double, long int, int, int, int, double, int, bool) except +
        void predefined(PyObject*, PyObject*, PyObject*, double, double, int,
                        double, double, double, long int, int, int, int, double, int) except +
