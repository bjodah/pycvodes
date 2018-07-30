# -*- mode: cython -*-
# -*- coding: utf-8 -*-

from anyode_numpy cimport NPY_TYPES
from cpython.ref cimport PyObject
from cvodes_cxx cimport realtype, indextype
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "anyode/anyode.hpp" namespace "AnyODE":
    cdef cppclass Info:
        unordered_map[string, int] nfo_int
        unordered_map[string, double] nfo_dbl
        unordered_map[string, vector[double]] nfo_vecdbl
        unordered_map[string, vector[int]] nfo_vecint

cdef extern from "cvodes_anyode_numpy.hpp" namespace "AnyODE":

    cdef cppclass CvodesPyOdeSys:
        CvodesPyOdeSys(indextype, PyObject*, PyObject*, PyObject*, PyObject*, PyObject*, PyObject*, int, int, int, int,
                       PyObject*, PyObject*, indextype)
        indextype get_ny()
        indextype get_nnz()
        NPY_TYPES get_int_type_tag()
        int get_nquads()
        int get_nroots()
        realtype get_dx0(realtype, realtype *) except +
        realtype get_dx_max(realtype, realtype *) except +
        bool autonomous_exprs
        bool use_get_dx_max
        bool record_rhs_xvals
        bool record_jac_xvals
        bool record_order
        bool record_fpe
        bool record_steps
        int mlower, mupper, nroots
        Info current_info
        int nfev, njev, njvev
        indextype nnz
        void * integrator