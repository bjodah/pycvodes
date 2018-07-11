# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType, LinSol

cdef extern from "cvodes_anyode_parallel.hpp" namespace "cvodes_anyode_parallel":
    cdef vector[pair[int, vector[int]]] multi_adaptive[U](
        double **,
        int *,
        vector[U*],
        vector[double],
        double,
        LMM,
        const double * const,
        long int,
        const double *,
        const double *,
        const double *,
        bool,
        IterType,
        LinSol,
        int,
        double,
        unsigned,
        bool,
        int,
        bool,
        bool
    ) except +

    cdef vector[pair[int, pair[vector[int], vector[double]]]] multi_predefined[U](
        vector[U*],
        vector[double],
        double,
        LMM,
        const double * const,
        size_t,
        const double * const,
        double * const,
        long int,
        double *,
        double *,
        double *,
        bool,
        IterType,
        LinSol,
        int,
        double,
        unsigned,
        int,
        bool,
        bool
    ) except +
