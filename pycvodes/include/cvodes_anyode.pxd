# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType

cdef extern from "cvodes_anyode.hpp" namespace "cvodes_anyode":
    cdef int simple_adaptive[U](
        double **,
        int *,
        U * const,
        vector[double],
        double,
        LMM,
        const double,
        vector[int]&,
        long int,
        double,
        double,
        double,
        bool,
        IterType,
        int,
        int,
        double,
        unsigned,
        bool,
        int,
        bool,
        bool,
        int,
        double **
    ) except +

    cdef int simple_predefined[U](
        U * const,
        vector[double],
        double,
        LMM,
        const double * const,
        size_t,
        const double * const,
        double * const,
        vector[int]&,
        vector[double]&,
        long int,
        double,
        double,
        double,
        bool,
        IterType,
        int,
        int,
        double,
        unsigned,
        int,
        bool,
        bool,
        double *
    ) except +
