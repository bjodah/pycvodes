# -*- coding: utf-8; mode: cython -*-

from cpython.ref cimport PyObject
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "anyode/anyode_numpy.hpp" namespace "AnyODE":
    cdef cppclass PyOdeSys:
        PyOdeSys(int, PyObject*, PyObject*, PyObject*, PyObject*, int, int, int, PyObject*, PyObject*)
        int get_ny()
        double get_dx0(double, double *) except +
        double get_dx_max(double, double *) except +
        bool use_get_dx_max
        int mlower, mupper, nroots
        unordered_map[string, int] last_integration_info
        unordered_map[string, double] last_integration_info_dbl
        int nfev, njev
        void * integrator
