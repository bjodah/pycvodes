# -*- coding: utf-8; mode: cython -*-

from cpython.ref cimport PyObject
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string


cdef extern from "anyode/anyode_numpy.hpp" namespace "AnyODE":
     cdef cppclass PyOdeSys:
         PyOdeSys(int, PyObject*, PyObject*, PyObject*, PyObject*, int, int, int)
         int get_ny()
         int mlower, mupper, nroots
         unordered_map[string, int] last_integration_info
         unordered_map[string, double] last_integration_info_dbl
         int nfev, njev
         void * integrator
