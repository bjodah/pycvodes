# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cvodes_cxx cimport LMM, IterType

cdef extern from "cvodes_anyode.hpp" namespace "cvodes_anyode":
    cdef pair[vector[double], vector[double]] simple_adaptive[U](
        U * const,
        vector[double],
        double,
        LMM,
        const double * const,
        double,
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
        bool
    ) except +

    cdef void simple_predefined[U](
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
        unsigned
    ) except +
