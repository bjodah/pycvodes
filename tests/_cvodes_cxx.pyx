from libcpp.vector cimport vector
from cvodes_cxx cimport simple_adaptive


cdef extern from "testing_utils.hpp":
    cppclass Decay:
        Decay(double)


cdef class PyDecay:
    cdef Decay *thisptr

    def __cinit__(self, double k):
        self.thisptr = new Decay(k)

    def __dealloc__(self):
        del self.thisptr

    def adaptive(self, double y0, double t):
        cdef:
            vector[int] root_indices
            int ADAMS = 1  # grep "define CV_ADAMS" cvodes.h
        return simple_adaptive[Decay](self.thisptr, [1e-10], 1e-10, ADAMS, &y0, 0.0, t, root_indices)
