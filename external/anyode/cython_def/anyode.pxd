# -*- coding: utf-8; mode: cython -*-

from cpython.ref cimport PyObject

cdef extern from "anyode/anyode.hpp" namespace "AnyODE":
     cdef cppclass OdeSysBase:
         int nfev, njev
