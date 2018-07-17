#pragma once

#include <Python.h>
#include <numpy/arrayobject.h>
#include <anyode/anyode_iterative.hpp>
#include <anyode/anyode_matrix.hpp> // DenseMatrix
#include <anyode/anyode_decomposition.hpp>  // DenseLU


BEGIN_NAMESPACE(AnyODE)
struct PyOdeSys : public AnyODE::OdeSysIterativeBase<double, int, DenseMatrix<double>, DenseLU<double>> {
    int ny;
    PyObject *py_rhs, *py_jac, *py_jtimes, *py_quads, *py_roots, *py_kwargs, *py_dx0cb, *py_dx_max_cb;
    int mlower, mupper, nquads, nroots;
    int nnz;
    PyOdeSys(int ny, PyObject * py_rhs, PyObject * py_jac=nullptr, PyObject * py_jtimes=nullptr,
             PyObject * py_quads=nullptr,
             PyObject * py_roots=nullptr, PyObject * py_kwargs=nullptr, int mlower=-1,
             int mupper=-1, int nquads=0, int nroots=0, PyObject * py_dx0cb=nullptr,
             PyObject * py_dx_max_cb=nullptr, int nnz=-1) :
        ny(ny), py_rhs(py_rhs), py_jac(py_jac), py_jtimes(py_jtimes),
        py_quads(py_quads), py_roots(py_roots),
        py_kwargs(py_kwargs), py_dx0cb(py_dx0cb), py_dx_max_cb(py_dx_max_cb),
        mlower(mlower), mupper(mupper), nquads(nquads), nroots(nroots),
        nnz(nnz)
    {
        if (py_rhs == nullptr){
            throw std::runtime_error("py_rhs must not be nullptr");
        }
        if ((py_dx_max_cb != nullptr) && (py_dx_max_cb != Py_None)) {
            this->use_get_dx_max = true;
        }
        Py_INCREF(py_rhs);
        Py_XINCREF(py_jac);
        Py_XINCREF(py_jtimes);
        Py_XINCREF(py_quads);
        Py_XINCREF(py_roots);
        if (py_kwargs == Py_None){
            Py_DECREF(Py_None);
            this->py_kwargs = nullptr;
        } else {
            Py_XINCREF(py_kwargs);
        }
    }
    virtual ~PyOdeSys() {
        Py_DECREF(py_rhs);
        Py_XDECREF(py_jac);
        Py_XDECREF(py_jtimes);
        Py_XDECREF(py_quads);
        Py_XDECREF(py_roots);
        Py_XDECREF(py_kwargs);
    }
    int get_ny() const override { return ny; }
    int get_mlower() const override { return mlower; }
    int get_mupper() const override { return mupper; }
    int get_nnz() const override { return nnz; }
    int get_nquads() const override { return nquads; }
    int get_nroots() const override { return nroots; }
    double get_dx0(double t, const double * const y) override {
        if (py_dx0cb == nullptr || py_dx0cb == Py_None) {
            return default_dx0;
        }
        npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(
            1, dims, type_tag, static_cast<void*>(const_cast<double*>(y)));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_arglist = Py_BuildValue("(dO)", (double)(t), py_yarr);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_dx0cb, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_yarr);
        if (py_result == nullptr) {
            throw std::runtime_error("get_dx0 failed (dx0cb failed)");
        }
        double res = PyFloat_AsDouble(py_result);
        Py_DECREF(py_result);
        if ((PyErr_Occurred()) && (res == -1.0)) {
            throw std::runtime_error("get_dx0 failed (value returned by dx0cb could not be converted to float)");
        }
        return res;
    }
    double get_dx_max(double t, const double * const y) override {
        if (py_dx_max_cb == nullptr || py_dx_max_cb == Py_None) {
            return INFINITY;
        }
        npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(
            1, dims, type_tag, static_cast<void*>(const_cast<double*>(y)));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_arglist = Py_BuildValue("(dO)", (double)(t), py_yarr);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_dx_max_cb, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_yarr);
        if (py_result == nullptr) {
            throw std::runtime_error("get_dx_max failed (dx_max_cb failed)");
        }
        double res = PyFloat_AsDouble(py_result);
        Py_DECREF(py_result);
        if (PyErr_Occurred() && (res == -1.0)) {
            throw std::runtime_error("get_dx_max failed (value returned by dx_max_cb could not be converted to float)");
        }
        return res;
    }
    Status handle_status_(PyObject * py_result, const std::string what_arg){
        if (py_result == nullptr){
            throw std::runtime_error(what_arg + " failed");
        } else if (py_result == Py_None){
            Py_DECREF(py_result);
            return AnyODE::Status::success;
        }
        long result = PyInt_AsLong(py_result);
        Py_DECREF(py_result);


        if ((PyErr_Occurred() && (result == -1)) ||
            (result == static_cast<long int>(AnyODE::Status::unrecoverable_error))) {
            return AnyODE::Status::unrecoverable_error;
        } else if (result == static_cast<long int>(AnyODE::Status::recoverable_error)) {
            return AnyODE::Status::recoverable_error;
        } else if (result == static_cast<long int>(AnyODE::Status::success)) {
            return AnyODE::Status::success;
        }
        throw std::runtime_error(what_arg + " did not return None, -1, 0 or 1");
    }
    Status rhs(double t, const double * const y, double * const dydt) override {
        npy_intp dims[1] { static_cast<npy_intp>(this->ny) } ;
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(
            1, dims, type_tag, static_cast<void*>(const_cast<double*>(y)));
        PyObject * py_dydt = PyArray_SimpleNewFromData(
            1, dims, type_tag, static_cast<void*>(dydt));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_arglist = Py_BuildValue("(dOO)", (double)(t), py_yarr, py_dydt);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_rhs, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_dydt);
        Py_DECREF(py_yarr);
        this->nfev++;
        return handle_status_(py_result, "rhs");
    }
    AnyODE::Status quads(double t, const double * const y, double * const out) override {
        npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
        npy_intp rdims[1] { static_cast<npy_intp>(this->get_nquads()) };
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(
            1, ydims, type_tag, static_cast<void*>(const_cast<double*>(y)));
        PyObject * py_out = PyArray_SimpleNewFromData(
            1, rdims, type_tag, static_cast<void*>(out));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_arglist = Py_BuildValue("(dOO)", t, py_yarr, py_out);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_quads, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_out);
        Py_DECREF(py_yarr);
        return handle_status_(py_result, "quads");
    }
    AnyODE::Status roots(double t, const double * const y, double * const out) override {
        npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
        npy_intp rdims[1] { static_cast<npy_intp>(this->get_nroots()) };
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(
            1, ydims, type_tag, static_cast<void*>(const_cast<double*>(y)));
        PyObject * py_out = PyArray_SimpleNewFromData(
            1, rdims, type_tag, static_cast<void*>(out));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_arglist = Py_BuildValue("(dOO)", t, py_yarr, py_out);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_roots, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_out);
        Py_DECREF(py_yarr);
        return handle_status_(py_result, "roots");
    }
    AnyODE::Status call_py_jac(double t, const double * const y, const double * const fy,
                               PyObject * py_jmat, double * const dfdt){
        npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *>(y));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_dfdt = (dfdt == nullptr) ? Py_BuildValue("") : PyArray_SimpleNewFromData(
            1, ydims, type_tag, static_cast<void*>(dfdt));
        PyObject * py_fy;
        if (fy) {
            py_fy = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *>(fy));
            PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_fy), NPY_ARRAY_WRITEABLE);  // make fy read-only
        } else {
            py_fy = Py_BuildValue(""); // Py_None with incref
        }
        // Call jac with signature: (t, y[:], Jmat[:, :], dfdt[:]=None, fy[:]=None)
        // (NumPy takes cares of row vs. column major ordering. User responsible for dense/banded.)
        PyObject * py_arglist = Py_BuildValue("(dOOOO)", (double)t, py_yarr, py_jmat, py_dfdt, py_fy);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_jac, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_fy);
        Py_DECREF(py_dfdt);
        Py_DECREF(py_yarr);
        this->njev++;
        return handle_status_(py_result, "jac");
    }
    AnyODE::Status jtimes(const double * const v, double * const Jv,
                          double x, const double * const y, const double * const fy) override {
        npy_intp ydims[1] { static_cast<npy_intp>(this->ny) };
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_yarr = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *>(y));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_varr = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *>(v));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_varr), NPY_ARRAY_WRITEABLE);  // make varr read-only
        PyObject * py_Jv = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *> (Jv));
        PyObject * py_fy;
        if (fy) {
            py_fy = PyArray_SimpleNewFromData(1, ydims, type_tag, const_cast<double *>(fy));
            PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_fy), NPY_ARRAY_WRITEABLE);  // make fy read-only
        } else {
            py_fy = Py_BuildValue(""); // Py_None with incref
        }
        // Call jtimes with signature: (v[:], Jv[:], x, y[:], fy[:])
        PyObject * py_arglist = Py_BuildValue("(OOdOO)", py_varr, py_Jv, (double) x, py_yarr, py_fy);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_jtimes, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_Jv);
        Py_DECREF(py_fy);
        Py_DECREF(py_yarr);
        Py_DECREF(py_varr);
        this->njvev++;
        return handle_status_(py_result, "jtimes");
    }
    AnyODE::Status dense_jac_cmaj(double t, const double * const y, const double * const fy,
                                  double * const jac, long int ldim, double * const dfdt=nullptr) override {
        npy_intp Jdims[2] { static_cast<npy_intp>(this->ny), static_cast<npy_intp>(this->ny) };
        npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
        int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE;
        if (ldim == Jdims[0]) {
            flags |= NPY_ARRAY_F_CONTIGUOUS;
        }
        const auto type_tag = NPY_DOUBLE;
        PyObject * py_jmat = PyArray_New(
            &PyArray_Type, 2, Jdims, type_tag, strides,
            static_cast<void *>(const_cast<double *>(jac)), sizeof(double),
            flags, nullptr);
        AnyODE::Status status = call_py_jac(t, y, fy, py_jmat, dfdt);
        Py_DECREF(py_jmat);
        return status;
    }
    AnyODE::Status dense_jac_rmaj(double t, const double * const y, const double * const fy,
                                  double * const jac, long int ldim, double * const dfdt=nullptr) override {
        npy_intp Jdims[2] { static_cast<npy_intp>(this->ny), static_cast<npy_intp>(this->ny) };
        npy_intp strides[2] { static_cast<npy_intp>(ldim*sizeof(double)), sizeof(double) };
        const auto type_tag = NPY_DOUBLE;
        int flags = NPY_ARRAY_ALIGNED| NPY_ARRAY_WRITEABLE;
        if (ldim == Jdims[1]) {
            flags |= NPY_ARRAY_C_CONTIGUOUS;
        }
        PyObject * py_jmat = PyArray_New(
            &PyArray_Type, 2, Jdims, type_tag, strides,
            static_cast<void *>(const_cast<double *>(jac)), sizeof(double), flags, nullptr);
        AnyODE::Status status = call_py_jac(t, y, fy, py_jmat, dfdt);
        Py_DECREF(py_jmat);
        return status;
    }
    AnyODE::Status sparse_jac_csc(double t, const double * const y, const double * const fy,
                                  double * const data, int * const colptrs, int * const rowvals) override {
        npy_intp y_dims[1] { static_cast<npy_intp>(this->ny) };
        npy_intp data_dims[1] { static_cast<npy_intp>(this->nnz) };
        npy_intp colptrs_dims[1] { static_cast<npy_intp>(this->ny + 1) };

        PyObject * py_yarr = PyArray_SimpleNewFromData(1, y_dims, NPY_DOUBLE, const_cast<double *>(y));
        PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_yarr), NPY_ARRAY_WRITEABLE);  // make yarr read-only
        PyObject * py_fy;
        if (fy) {
            py_fy = PyArray_SimpleNewFromData(1, y_dims, NPY_DOUBLE, const_cast<double *>(fy));
            PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(py_fy), NPY_ARRAY_WRITEABLE);  // make fy read-only
        } else {
            py_fy = Py_BuildValue(""); // Py_None with incref
        }
        PyObject * py_data = PyArray_SimpleNewFromData(1, data_dims, NPY_DOUBLE, static_cast<double *>(data));
        PyObject * py_colptrs = PyArray_SimpleNewFromData(1, colptrs_dims, NPY_INT, static_cast<int *>(colptrs));
        PyObject * py_rowvals = PyArray_SimpleNewFromData(1, data_dims, NPY_INT, static_cast<int *>(rowvals));

        // Call sparse jac with signature: (t, y[:], data[:], colptrs[:], rowvals[:]
        PyObject * py_arglist = Py_BuildValue("(dOOOO)", (double) t, py_yarr, py_data, py_colptrs, py_rowvals);
        PyObject * py_result = PyEval_CallObjectWithKeywords(this->py_jac, py_arglist, this->py_kwargs);
        Py_DECREF(py_arglist);
        Py_DECREF(py_fy);
        Py_DECREF(py_yarr);
        Py_DECREF(py_data);
        Py_DECREF(py_colptrs);
        Py_DECREF(py_rowvals);
        this->njev++;
        return handle_status_(py_result, "jac");
    }
    AnyODE::Status banded_jac_cmaj(double t, const double * const y, const double * const fy,
                                   double * const jac, long int ldim) override {
        npy_intp Jdims[2] { 1 + this->mlower + this->mupper, static_cast<npy_intp>(this->ny) };
        npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
        const auto type_tag = NPY_DOUBLE;
        int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE;
        if (ldim == Jdims[0] ) {
            flags |= NPY_ARRAY_F_CONTIGUOUS;
        }
        PyObject * py_jmat = PyArray_New(
            &PyArray_Type, 2, Jdims, type_tag, strides,
            static_cast<void *>(const_cast<double *>(jac)), sizeof(double), flags, nullptr);
        AnyODE::Status status = call_py_jac(t, y, fy, py_jmat, nullptr);
        Py_DECREF(py_jmat);
        return status;
    }
};
END_NAMESPACE(AnyODE)
