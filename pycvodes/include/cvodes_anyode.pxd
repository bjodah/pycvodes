# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType, LinSol, realtype

cdef extern from "cvodes_anyode.hpp" namespace "cvodes_anyode":
    cdef int simple_adaptive[U](
        realtype **,
        int *,
        U * const,
        vector[realtype],
        const realtype,
        const LMM,
        const realtype,
        vector[int]&,
        long int,
        const realtype,
        const realtype,
        const realtype,
        const bool,
        IterType,
        LinSol,
        const int,
        const realtype,
        const unsigned,
        const bool,
        const int,
        const bool,
        const bool,
        int,
        realtype **,
        vector[double]&
    ) except +

    cdef int simple_predefined[U](
        U * const,
        vector[realtype],
        const realtype,
        const LMM,
        const realtype * const,
        const size_t,
        const realtype * const,
        realtype * const,
        vector[int]&,
        vector[realtype]&,
        const long int,
        realtype,
        const realtype,
        const realtype,
        const bool,
        IterType,
        LinSol,
        const int,
        const realtype,
        const unsigned,
        const int,
        const bool,
        const bool,
        realtype *,
	vector[double]&
    ) except +
