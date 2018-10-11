# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType, LinSol, realtype

cdef extern from "cvodes_anyode_parallel.hpp" namespace "cvodes_anyode_parallel":
    cdef vector[pair[int, vector[int]]] multi_adaptive[U](
        realtype **,
        int *,
        vector[U*],
        vector[realtype],
        const realtype,
        const LMM,
        const realtype *,
        const long int,
        const realtype *,
        const realtype *,
        const realtype *,
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
        realtype***,
        vector[realtype]&
    ) except +

    cdef vector[pair[int, pair[vector[int], vector[realtype]]]] multi_predefined[U](
        vector[U*],
        vector[realtype],
        const realtype,
        const LMM,
        realtype *,
        const size_t,
        realtype *,
        realtype *,
        const long int,
        const realtype *,
        const realtype *,
        const realtype *,
        const bool,
        IterType,
        LinSol,
        const int,
        const realtype,
        const unsigned,
        const int,
        const bool,
        const bool,
        realtype **,
        vector[realtype]&
    ) except +
