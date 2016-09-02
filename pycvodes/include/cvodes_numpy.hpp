#pragma once
#include <Python.h>
#include <numpy/arrayobject.h>
#include <chrono>

#include <utility> // std::pair
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <string> // std::string

#include "anyode.hpp"
#include "cvodes_cxx.hpp"
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */


namespace cvodes_numpy{
    // typedef int (*RhsFn)(double t, const double y[], double dydt[], void *params);
    // typedef int (*JacFn)(double t, const double y[], double *dfdy, double dfdt[], void *params);
    // int rhs(double t, const double y[], double f[], void * params);
    // int jac(double x, const double y[], double *dfdy, double dfdt[], void *params);
    // int roots(double x, const double y[], double * const out, void * params);


    struct PyCvodes : public AnyODE::OdeSysBase {
        PyObject *py_rhs, *py_jac, *py_roots;
        const size_t ny;
        size_t nfev, njev;
        double time_cpu, time_wall;
        const int mlower, mupper, nroots;
        std::vector<realtype> xout;
        std::vector<realtype> yout;
        std::vector<int> root_indices;
        std::vector<double> roots_output;

        PyCvodes(PyObject * py_rhs, PyObject * py_jac, PyObject * py_roots, size_t ny, int ml=-1, int mu=-1, int nroots=0) :
            py_rhs(py_rhs), py_jac(py_jac), py_roots(py_roots),
            ny(ny), mlower(ml), mupper(mu), nroots(nroots) {}

        int get_ny() const override { return this->ny; }
        int get_mlower() const override { return this->mlower; }
        int get_mupper() const override { return this->mupper; }
        int get_nroots() const override { return this->nroots; }

        size_t adaptive(PyObject *py_y0, realtype x0, realtype xend,
                        realtype atol, realtype rtol, int step_type_idx,
                        realtype dx0, realtype dx_min=0.0, realtype dx_max=0.0, long int mxsteps=0,
                        int iter_type=0, int linear_solver=0, int maxl=5, realtype eps_lin=0.0, int nderiv=0,
                        bool return_on_root=false){
            std::clock_t cputime0 = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            const bool with_jacobian = py_jac != Py_None;
            auto y0 = (realtype*)PyArray_GETPTR1(py_y0, 0);
            nfev = 0; njev = 0;
            this->root_indices.clear();
            auto xy_out = cvodes_cxx::simple_adaptive(this, std::vector<realtype>(1, atol), rtol, step_type_idx, y0, x0,
                                                      xend, this->root_indices, dx0, dx_min, dx_max, mxsteps, with_jacobian,
                                                      iter_type, linear_solver, maxl, eps_lin, nderiv, return_on_root);
            this->xout = xy_out.first;
            this->yout = xy_out.second;
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            this->time_wall = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();
            return this->xout.size();
        }

        void predefined(PyObject *py_y0, PyObject *py_xout, PyObject *py_yout,
                        realtype atol, realtype rtol,
                        int step_type_idx,
                        realtype dx0, realtype dx_min=0.0, realtype dx_max=0.0, long int mxsteps=0,
                        int iter_type=0, int linear_solver=0, int maxl=5, realtype eps_lin=0.0, int nderiv=0) {
            std::clock_t cputime0 = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();
            auto y0 = (realtype*)PyArray_GETPTR1(py_y0, 0);
            auto xout = (realtype*)PyArray_GETPTR1(py_xout, 0);
            auto yout = (realtype*)PyArray_GETPTR1(py_yout, 0);
            const npy_intp nt = PyArray_DIMS(py_xout)[0];
            const bool with_jacobian = py_jac != Py_None;
            nfev = 0; njev = 0;
            this->root_indices.clear();
            this->roots_output.clear();
            cvodes_cxx::simple_predefined<PyCvodes>(this, std::vector<realtype>(1, atol), rtol,
                                                    step_type_idx, y0, nt, xout, yout, this->root_indices,
                                                    this->roots_output, dx0, dx_min,
                                                    dx_max, mxsteps, with_jacobian, iter_type, linear_solver, maxl,
                                                    eps_lin, nderiv);
            this->time_cpu = (std::clock() - cputime0) / (double)CLOCKS_PER_SEC;
            this->time_wall = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();
        }
        AnyODE::Status handle_status_(PyObject * py_result, const std::string what_arg){
            if (py_result == nullptr){
                throw std::runtime_error(what_arg + " failed");
            } else if (py_result == Py_None){
                Py_DECREF(py_result);
                return AnyODE::Status::success;
            }
            long result = PyInt_AsLong(py_result);
            Py_DECREF(py_result);
            if ((PyErr_Occurred() and result == -1) or result == static_cast<long int>(AnyODE::Status::unrecoverable_error))
                return AnyODE::Status::unrecoverable_error;
            else if (result == static_cast<long int>(AnyODE::Status::recoverable_error))
                return AnyODE::Status::recoverable_error;
            else if (result == static_cast<long int>(AnyODE::Status::success))
                return AnyODE::Status::success;
            throw std::runtime_error(what_arg + " did not return None, -1, 0 or 1");
        }
        virtual AnyODE::Status rhs(realtype xval, const realtype * const y, realtype * const dydx){
            npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
            const auto type_tag = (sizeof(realtype) == 8) ? NPY_DOUBLE : NPY_LONGDOUBLE;
            PyObject * py_yarr = PyArray_SimpleNewFromData(1, dims, type_tag, static_cast<void*>(const_cast<realtype*>(y)));
            PyObject * py_dydx = PyArray_SimpleNewFromData(1, dims, type_tag, static_cast<void*>(dydx));
            PyObject * py_arglist = Py_BuildValue("(dOO)", (double)(xval), py_yarr, py_dydx);
            PyObject * py_result = PyEval_CallObject(this->py_rhs, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_dydx);
            Py_DECREF(py_yarr);
            nfev++;
            return handle_status_(py_result, "rhs");
        }
        virtual AnyODE::Status roots(realtype xval, const realtype * const y, realtype * const out){
            npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
            npy_intp rdims[1] { static_cast<npy_intp>(this->get_nroots()) };
            const auto type_tag = (sizeof(realtype) == 8) ? NPY_DOUBLE : NPY_LONGDOUBLE;
            PyObject * py_yarr = PyArray_SimpleNewFromData(
                1, ydims, type_tag, static_cast<void*>(const_cast<realtype*>(y)));
            PyObject * py_out = PyArray_SimpleNewFromData(
                1, rdims, type_tag, static_cast<void*>(out));
            PyObject * py_arglist = Py_BuildValue("(dOO)", xval, py_yarr, py_out);
            PyObject * py_result = PyEval_CallObject(this->py_roots, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_out);
            Py_DECREF(py_yarr);
            return handle_status_(py_result, "roots");
        }
        AnyODE::Status call_py_jac(realtype t, const realtype * const y, const realtype * const fy,
                         PyObject * py_jmat){
            npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
            const auto type_tag = (sizeof(realtype) == 8) ? NPY_DOUBLE : NPY_LONGDOUBLE;
            PyObject * py_yarr = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<realtype *>(y));
            PyObject * py_fy = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<realtype *>(fy));
            PyObject * py_arglist = Py_BuildValue("(dOOOO)", (double)t, py_yarr, py_jmat, Py_None,py_fy);
            PyObject * py_result = PyEval_CallObject(this->py_jac, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_fy);
            Py_DECREF(py_yarr);
            njev++;
            return handle_status_(py_result, "jac");
        }
        virtual AnyODE::Status dense_jac_cmaj(realtype t, const realtype * const y, const realtype * const fy,
                                      realtype * const jac, long int ldim){
            npy_intp Jdims[2] { static_cast<npy_intp>(this->ny), static_cast<npy_intp>(this->ny) };
            npy_intp strides[2] { sizeof(realtype), static_cast<npy_intp>(ldim*sizeof(realtype)) };
            const auto type_tag = (sizeof(realtype) == 8) ? NPY_DOUBLE : NPY_LONGDOUBLE;
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, type_tag, strides,
                static_cast<void *>(const_cast<realtype *>(jac)), sizeof(realtype),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            AnyODE::Status status = call_py_jac(t, y, fy, py_jmat);
            Py_DECREF(py_jmat);
            return status;
        }
        virtual AnyODE::Status banded_jac_cmaj(realtype t, const realtype * const y, const realtype * const fy,
                                    realtype * const jac, long int ldim){
            npy_intp Jdims[2] { 1 + this->mlower + this->mupper, static_cast<npy_intp>(this->ny) };
            npy_intp strides[2] { sizeof(realtype), static_cast<npy_intp>(ldim*sizeof(realtype)) };
            const auto type_tag = (sizeof(realtype) == 8) ? NPY_DOUBLE : NPY_LONGDOUBLE;
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, type_tag, strides,
                static_cast<void *>(const_cast<realtype *>(jac)), sizeof(realtype),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            AnyODE::Status status = call_py_jac(t, y, fy, py_jmat);
            Py_DECREF(py_jmat);
            return status;
        }
        virtual AnyODE::Status jac_times_vec(const realtype * const vec,
                                             realtype * const out,
                                             realtype t,
                                             const realtype * const y,
                                             const realtype * const fy
                                             )
        {
            cvodes_cxx::ignore(vec);
            cvodes_cxx::ignore(out);
            cvodes_cxx::ignore(t);
            cvodes_cxx::ignore(y);
            cvodes_cxx::ignore(fy);
            throw std::runtime_error("Not implemented!");
        }
        virtual AnyODE::Status prec_setup(realtype t,
                                          const realtype * const __restrict__ y,
                                          const realtype * const __restrict__ fy,
                                          bool jok,
                                          bool& jac_recomputed,
                                          realtype gamma)
        {
            cvodes_cxx::ignore(t);
            cvodes_cxx::ignore(y);
            cvodes_cxx::ignore(fy);
            cvodes_cxx::ignore(jok);
            cvodes_cxx::ignore(jac_recomputed);
            cvodes_cxx::ignore(gamma);
            throw std::runtime_error("Not implemented!");
        }
        virtual AnyODE::Status prec_solve_left(const realtype t,
                                               const realtype * const __restrict__ y,
                                               const realtype * const __restrict__ fy,
                                               const realtype * const __restrict__ r,
                                               realtype * const __restrict__ z,
                                               realtype gamma,
                                               realtype delta,
                                               const realtype * const __restrict__ ewt)
        {
            cvodes_cxx::ignore(t);
            cvodes_cxx::ignore(y);
            cvodes_cxx::ignore(fy);
            cvodes_cxx::ignore(r);
            cvodes_cxx::ignore(z);
            cvodes_cxx::ignore(gamma);
            cvodes_cxx::ignore(delta);
            cvodes_cxx::ignore(ewt);
            throw std::runtime_error("Not implemented!");
        }
    };
}
