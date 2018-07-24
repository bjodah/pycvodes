# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from cvodes_cxx import realtype, indextype
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType, LinSol

cdef extern from "cvodes_anyode_parallel.hpp" namespace "cvodes_anyode_parallel":
    cdef vector[pair[int, vector[int]]] multi_adaptive[U](
        realtype **,
        int *,
        vector[U*],
        vector[realtype],
        realtype,
        LMM,
        const realtype *,
        long int,
        const realtype *,
        const realtype *,
        const realtype *,
        bool,
        IterType,
        LinSol,
        int,
        realtype,
        unsigned,
        bool,
        int,
        bool,
        bool
    ) except +

    cdef vector[pair[int, pair[vector[int], vector[realtype]]]] multi_predefined[U](
        vector[U*],
        const vector[realtype],
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
        bool,
        IterType,
        LinSol,
        const int,
        const realtype,
        const unsigned,
        int,
        bool,
        bool
    ) except +
