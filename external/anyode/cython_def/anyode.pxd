# -*- coding: utf-8; mode: cython -*-

cdef extern from "anyode/anyode.hpp" namespace "AnyODE":
     cdef cppclass OdeSysBase:
         int nfev, njev
