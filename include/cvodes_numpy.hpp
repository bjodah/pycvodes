#ifndef CVODES_NUMPY_HPP_MRTG5EIRAJFXZBUE7MHMKQUGYA
#define CVODES_NUMPY_HPP_MRTG5EIRAJFXZBUE7MHMKQUGYA
#include <Python.h>
#include <numpy/arrayobject.h>
#include <chrono>

#include <utility> // std::pair
#include <vector> // std::vector

#include "cvodes_cxx.hpp"
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */


namespace cvodes_numpy{
    // typedef int (*RhsFn)(double t, const double y[], double dydt[], void *params);
    // typedef int (*JacFn)(double t, const double y[], double *dfdy, double dfdt[], void *params);
    // int rhs(double t, const double y[], double f[], void * params);
    // int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);
    // int roots(double x, const double y[], double * const out, void * params);


    class PyCvodes {
    public:
        PyObject *py_rhs, *py_jac, *py_roots;
        const size_t ny;
        size_t nrhs, njac;
        double time_cpu;
        const int mlower, mupper, nroots;
        std::vector<double> xout;
        std::vector<double> yout;
        std::vector<int> root_indices;

        PyCvodes(PyObject * py_rhs, PyObject * py_jac, PyObject * py_roots, size_t ny, int ml=-1, int mu=-1, int nroots=0) :
            py_rhs(py_rhs), py_jac(py_jac), py_roots(py_roots), ny(ny), mlower(ml), mupper(mu), nroots(nroots) {}

        size_t adaptive(PyObject *py_y0, double x0, double xend,
                        double atol, double rtol, int step_type_idx,
                        double dx0, double dx_min=0.0, double dx_max=0.0, long int mxsteps=0,
                        int iterative=0, int nderiv=0, int sparse=0, bool return_on_root=false){
            std::clock_t cputime0 = std::clock();
            const bool with_jacobian = py_jac != Py_None;
            auto y0 = (double*)PyArray_GETPTR1(py_y0, 0);
            nrhs = 0; njac = 0;
            this->root_indices.clear();
            auto xy_out = cvodes_cxx::simple_adaptive(this, std::vector<double>(1, atol),
                                                      rtol, step_type_idx, y0, x0, xend, dx0, dx_min, dx_max, mxsteps,
                                                      (this->mlower > -1) ? 2 : 1, with_jacobian,
                                                      iterative, nderiv, this->root_indices, sparse, return_on_root);
            this->xout = xy_out.first;
            this->yout = xy_out.second;
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            return this->xout.size();
        }

        void predefined(PyObject *py_y0, PyObject *py_xout, PyObject *py_yout,
                        double atol, double rtol,
                        int step_type_idx,
                        double dx0, double dx_min=0.0, double dx_max=0.0, long int mxsteps=0,
                        int iterative=0, int nderiv=0) {
            std::clock_t cputime0 = std::clock();
            auto y0 = (double*)PyArray_GETPTR1(py_y0, 0);
            auto xout = (double*)PyArray_GETPTR1(py_xout, 0);
            auto yout = (double*)PyArray_GETPTR1(py_yout, 0);
            const npy_intp nt = PyArray_DIMS(py_xout)[0];
            const bool with_jacobian = py_jac != Py_None;
            nrhs = 0; njac = 0;
            this->root_indices.clear();
            cvodes_cxx::simple_predefined<PyCvodes>(this, std::vector<double>(1, atol), rtol,
                                                    step_type_idx, y0, nt, xout, yout, dx0, dx_min, dx_max, mxsteps,
                                                    (this->mlower > -1) ? 2 : 1, with_jacobian,
                                                    iterative, nderiv, this->root_indices);
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
        }

        void rhs(double xval, const double * const y, double * const dydx){
            npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
            PyObject * py_yarr = PyArray_SimpleNewFromData(
                1, dims, NPY_DOUBLE, static_cast<void*>(const_cast<double*>(y)));
            PyObject * py_dydx = PyArray_SimpleNewFromData(
                1, dims, NPY_DOUBLE, static_cast<void*>(dydx));
            PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_dydx);
            PyObject * py_result = PyEval_CallObject(this->py_rhs, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_dydx);
            Py_DECREF(py_yarr);
            nrhs++;
            if (py_result == nullptr){
                throw std::runtime_error("rhs() failed");
            } else if (py_result != Py_None){
                // py_result is not None
                throw std::runtime_error("rhs() failed");
            }
            Py_DECREF(py_result);
        }
        void roots(double xval, const double * const y, double * const out){
            npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
            npy_intp rdims[1] { static_cast<npy_intp>(this->nroots) };
            PyObject * py_yarr = PyArray_SimpleNewFromData(
                1, ydims, NPY_DOUBLE, static_cast<void*>(const_cast<double*>(y)));
            PyObject * py_out = PyArray_SimpleNewFromData(
                1, rdims, NPY_DOUBLE, static_cast<void*>(out));
            PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_out);
            PyObject * py_result = PyEval_CallObject(this->py_roots, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_out);
            Py_DECREF(py_yarr);
            if (py_result == nullptr){
                throw std::runtime_error("roots() failed");
            } else if (py_result != Py_None){
                // py_result is not None
                throw std::runtime_error("roots() did not return None");
            }
            Py_DECREF(py_result);
        }
        void call_py_jac(double t, const double * const y, const double * const fy,
                         PyObject * py_jmat){
            npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
            PyObject * py_yarr = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE,
                const_cast<double *>(y));
            PyObject * py_fy = PyArray_SimpleNewFromData(1, ydims, NPY_DOUBLE,
                                                         const_cast<double *>(fy));
            PyObject * py_arglist = Py_BuildValue("(dOOOO)", t, py_yarr, py_jmat, Py_None,py_fy);
            PyObject * py_result = PyEval_CallObject(this->py_jac, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_fy);
            Py_DECREF(py_yarr);
            njac++;
            if (py_result == nullptr){
                throw std::runtime_error("jac() failed");
            } else if (py_result != Py_None){
                // py_result is not None
                throw std::runtime_error("jac() did not return None");
            }
            Py_DECREF(py_result);
        }
        void dense_jac_cmaj(double t, const double * const y, const double * const fy,
                            double * const jac, long int ldim){
            npy_intp Jdims[2] { static_cast<npy_intp>(this->ny), static_cast<npy_intp>(this->ny) };
            npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, NPY_DOUBLE, strides,
                static_cast<void *>(const_cast<double *>(jac)), sizeof(double),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            call_py_jac(t, y, fy, py_jmat);
            Py_DECREF(py_jmat);
        }
        void banded_padded_jac_cmaj(double t, const double * const y, const double * const fy,
                                    double * const jac, long int ldim){
            npy_intp Jdims[2] { 1 + this->mlower + this->mupper, static_cast<npy_intp>(this->ny) };
            npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, NPY_DOUBLE, strides,
                static_cast<void *>(const_cast<double *>(jac + this->mupper)), sizeof(double),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            call_py_jac(t, y, fy, py_jmat);
            Py_DECREF(py_jmat);
        }
        void jac_times_vec(const double * const vec,
                           double * const out,
                           double t,
                           const double * const y,
                           const double * const fy
                           )
        { throw std::runtime_error("Not implemented!"); }
        void prec_setup(double t,
                const double * const __restrict__ y,
                const double * const __restrict__ fy,
                bool jok, bool& jac_recomputed, double gamma)
        { throw std::runtime_error("Not implemented!"); }
        void prec_solve_left(const double t,
                        const double * const __restrict__ y,
                        const double * const __restrict__ fy,
                        const double * const __restrict__ r,
                        double * const __restrict__ z,
                        double gamma)
        { throw std::runtime_error("Not implemented!"); }
    };

} // namespace cvodes_numpy

#endif /* CVODES_NUMPY_HPP_MRTG5EIRAJFXZBUE7MHMKQUGYA */
