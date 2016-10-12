# -*- coding: utf-8; mode: cython -*-

from libcpp.string cimport string


cdef extern from "cvodes_cxx.hpp" namespace "cvodes_cxx":
    cdef cppclass LMM:
        pass  # LMM is an enum class

    cdef cppclass IterType:
        pass  # IterType is an enum class

    cdef cppclass CVodeIntegrator:
        CVodeIntegrator(LMM, IterType)

    cdef LMM lmm_from_name(string) nogil except +
    cdef IterType iter_type_from_name(string) nogil except +

cdef extern from "cvodes_cxx.hpp" namespace "cvodes_cxx::LMM":
    cdef LMM Adams
    cdef LMM BDF

cdef extern from "cvodes_cxx.hpp" namespace "cvodes_cxx::IterType":
    cdef IterType Functional
    cdef IterType Newton
