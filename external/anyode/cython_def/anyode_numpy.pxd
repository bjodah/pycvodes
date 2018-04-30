# -*- coding: utf-8; mode: cython -*-

from cpython.ref cimport PyObject
from libcpp cimport bool
from anyode cimport Info

cdef extern from "anyode/anyode_numpy.hpp" namespace "AnyODE":
    cdef cppclass PyOdeSys:
        PyOdeSys(int, PyObject*, PyObject*, PyObject*, PyObject*, PyObject*, int, int, int, int, PyObject*, PyObject*)
        int get_ny()
        int get_nquads()
        int get_nroots()
        double get_dx0(double, double *) except +
        double get_dx_max(double, double *) except +
        bool autonomous_exprs
        bool use_get_dx_max
        bool record_rhs_xvals
        bool record_jac_xvals
        bool record_order
        bool record_fpe
        bool record_steps
        int mlower, mupper, nroots
        Info current_info
        int nfev, njev
        void * integrator
