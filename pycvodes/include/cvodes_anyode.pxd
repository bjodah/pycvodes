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
        realtype,
        LMM,
        const realtype,
        vector[int]&,
        long int,
        realtype,
        realtype,
        realtype,
        bool,
        IterType,
        LinSol,
        int,
        realtype,
        unsigned,
        bool,
        int,
        bool,
        bool,
        int,
        realtype **
    ) except +

    cdef int simple_predefined[U](
        U * const,
        vector[realtype],
        realtype,
        LMM,
        const realtype * const,
        size_t,
        const realtype * const,
        realtype * const,
        vector[int]&,
        vector[realtype]&,
        long int,
        realtype,
        realtype,
        realtype,
        bool,
        IterType,
        LinSol,
        int,
        realtype,
        unsigned,
        int,
        bool,
        bool,
        realtype *
    ) except +
