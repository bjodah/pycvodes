# -*- coding: utf-8; mode: cython -*-

from libcpp cimport bool

cdef extern from "anyode/anyode.hpp" namespace "AnyODE":
     cdef cppclass OdeSysBase:
         int nfev, njev
         bool use_get_dx_max
