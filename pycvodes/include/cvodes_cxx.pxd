# -*- coding: utf-8 -*-
# -*- mode: cython -*-
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "cvodes_cxx.hpp" namespace "cvodes_cxx":
    cdef void simple_predefined[U](
        U * const, vector[double], double, int, const double * const, size_t, const double * const, double * const,
        vector[int]&, vector[double]&, double, double, double, long int, bool, int, int, int, double, int) except +

    cdef pair[vector[double], vector[double]] simple_adaptive[U](
        U * const, vector[double], double, int, const double * const, double, const double tend, vector[int]&, double,
        double, double, long int, bool, int, int, int, double, int, bool) except +
