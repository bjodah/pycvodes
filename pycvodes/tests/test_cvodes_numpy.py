# -*- coding: utf-8 -*-
from math import exp, pi
import os
import numpy as np
import pytest

from pycvodes import (
    integrate_adaptive, integrate_predefined, requires_jac, get_include
)


def test_get_include():
    assert get_include().endswith('include')
    assert 'cvodes_cxx.hpp' in os.listdir(get_include())


decay_analytic = {
    0: lambda y0, k, t: (
        y0[0] * np.exp(-k[0]*t)),
    1: lambda y0, k, t: (
        y0[1] * np.exp(-k[1] * t) + y0[0] * k[0] / (k[1] - k[0]) *
        (np.exp(-k[0]*t) - np.exp(-k[1]*t))),
    2: lambda y0, k, t: (
        y0[2] * np.exp(-k[2] * t) + y0[1] * k[1] / (k[2] - k[1]) *
        (np.exp(-k[1]*t) - np.exp(-k[2]*t)) +
        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
        (1 / (k[2] - k[0]) * (np.exp(-k[0]*t) - np.exp(-k[2]*t)) -
         1 / (k[2] - k[1]) * (np.exp(-k[1]*t) - np.exp(-k[2]*t))))
}


def decay_get_Cref(k, y0, tout):
    coeffs = list(k) + [0]*(3-len(k))
    return np.column_stack([
        decay_analytic[i](y0, coeffs, tout) for i in range(
            min(3, len(k)+1))])


def _get_f_j(k):
    k0, k1, k2 = k

    def f(t, y, fout):
        fout[0] = -k0*y[0]
        fout[1] = k0*y[0] - k1*y[1]
        fout[2] = k1*y[1] - k2*y[2]

    def j(t, y, jmat_out, dfdx_out=None, fy=None):
        jmat_out[0, 0] = -k0
        jmat_out[0, 1] = 0
        jmat_out[0, 2] = 0
        jmat_out[1, 0] = k0
        jmat_out[1, 1] = -k1
        jmat_out[1, 2] = 0
        jmat_out[2, 0] = 0
        jmat_out[2, 1] = k1
        jmat_out[2, 2] = -k2
        if dfdx_out is not None:
            dfdx_out[0] = 0
            dfdx_out[1] = 0
            dfdx_out[2] = 0
    return f, j

methods = [('adams', 1.8, False),
           ('adams', 1.8, True),
           ('bdf', 10, False),
           ('bdf', 9, True)]


def bandify(cb, mlower, mupper):
    def j(t, y, jmat_out, dfdx_out=None, fy=None):
        jmat = np.empty((len(y), len(y)))
        cb(t, y, jmat, dfdx_out, fy)
        for ci in range(len(y)):
            for ri in range(max(0, ci-mupper), min(len(y), ci+mlower+1)):
                jmat_out[ri - ci + mupper, ci] = jmat[ri, ci]
    return j


@pytest.mark.parametrize("method,forgiveness,banded", methods)
def test_integrate_predefined(method, forgiveness, banded):
    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j = _get_f_j(k)
    kwargs = {'method': method}
    if not use_jac:
        j = None
    else:
        if banded:
            j = bandify(j, 1, 0)
            kwargs['lband'] = 1
            kwargs['uband'] = 0
    xout = np.linspace(0, 3, 31)
    dx0 = 1e-10
    atol, rtol = 1e-8, 1e-8
    # Run twice to catch possible side-effects:
    yout, nfo = integrate_predefined(f, j, y0, xout, dx0, 1e-8, 1e-8, **kwargs)
    yout, nfo = integrate_predefined(f, j, y0, xout, dx0, 1e-8, 1e-8, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert nfo['nfev'] > 0
    assert nfo['time_cpu'] > 1e-9
    assert nfo['time_wall'] > 1e-9
    if method in requires_jac:
        assert nfo['njev'] > 0


@pytest.mark.parametrize("method,forgiveness,banded", methods)
def test_integrate_adaptive(method, forgiveness, banded):
    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(x0=0, xend=3, dx0=1e-10, atol=atol, rtol=rtol,
                  method=method, iter_type='newton')
    f, j = _get_f_j(k)
    if not use_jac:
        j = None
    else:
        if banded:
            j = bandify(j, 1, 0)
            kwargs['lband'] = 1
            kwargs['uband'] = 0
    # Run twice to catch possible side-effects:
    xout, yout, info = integrate_adaptive(f, j, y0, **kwargs)
    xout, yout, info = integrate_adaptive(f, j, y0, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=forgiveness*rtol,
                       atol=forgiveness*atol)
    assert info['nfev'] > 0
    if method in requires_jac:
        assert info['njev'] > 0


def test_derivative_1():
    def f(t, y, fout):
        fout[0] = y[0]
    kwargs = dict(dx0=0.0, atol=1e-8, rtol=1e-8, nderiv=1, method='adams')
    yout, info = integrate_predefined(f, None, [1], [0, 1, 2], **kwargs)
    assert yout.shape == (3, 2, 1)
    ref = np.array([
        [[exp(0)], [exp(0)]],
        [[exp(1)], [exp(1)]],
        [[exp(2)], [exp(2)]],
    ])
    assert np.allclose(yout, ref)


def test_derivative_2():
    def f(t, y, fout):
        fout[0] = y[0]
    kwargs = dict(dx0=0.0, atol=1e-12, rtol=1e-12, nderiv=3, method='adams')
    yout, info = integrate_predefined(f, None, [1], [0, 1, 2, 3, 4], **kwargs)
    assert yout.shape == (5, 4, 1)
    ref = np.array([
        [[exp(0)], [exp(0)], [0], [0]],  # higher order deriv. skipped for t0
        [[exp(1)]]*4,
        [[exp(2)]]*4,
        [[exp(3)]]*4,
        [[exp(4)]]*4,
    ])
    assert np.allclose(yout, ref)


def test_derivative_3():
    def f(t, y, fout):
        fout[0] = y[1]
        fout[1] = -y[0]
    kwargs = dict(dx0=0.0, atol=1e-13, rtol=1e-13, nderiv=2, method='adams',
                  iter_type='newton')
    xout, yout, info = integrate_adaptive(f, None, [0, 1], 0, 4*pi, **kwargs)
    assert yout.shape[1:] == (3, 2)
    sinx, cosx = np.sin(xout), np.cos(xout)
    ref = np.empty((len(xout), 3, 2))
    ref[:, 0, 0], ref[:, 0, 1] = sinx, cosx
    ref[:, 1, 0], ref[:, 1, 1] = cosx, -sinx
    ref[:, 2, 0], ref[:, 2, 1] = -sinx, -cosx
    discrepancy = yout[7:, ...] - ref[7:, ...]
    assert np.allclose(discrepancy, 0, rtol=1e-6, atol=1e-6)


def test_roots_adaptive():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    assert len(info['root_indices']) == 1
    assert np.min(np.abs(xout - 1)) < 1e-11


def test_roots_predefined():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1)
    xout = [0, 0.5, 1.5, 2]
    yout, info = integrate_predefined(f, None, [1], xout, **kwargs)
    discrepancy = np.exp(xout) - yout.flatten()
    assert np.allclose(discrepancy, 0)
    assert len(info['root_indices']) == 1
    assert info['root_indices'][0] == 2


def test_adaptive_nderiv_4():
    def f(t, y, fout):
        fout[0] = y[0]
    kwargs = dict(dx0=1e-4, atol=1e-4, rtol=1e-12, method='adams',
                  nderiv=4)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    discrepancy = np.exp(xout) - yout[:, 0].flatten()
    assert np.allclose(discrepancy, 0, atol=1e-3)


def test_adaptive_nderiv():
    def f(t, y, fout):
        fout[0] = y[0]
    kwargs = dict(dx0=1e-4, atol=1e-4, rtol=1e-12, method='adams',
                  nderiv=4)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    discrepancy = np.exp(xout) - yout[:, 0].flatten()
    assert np.allclose(discrepancy, 0, atol=1e-3)


def test_return_on_root():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1, return_on_root=True)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    assert len(info['root_indices']) == 1
    assert abs(xout[-1] - 1) < 1e-11
    assert abs(yout[-1, 0] - exp(1)) < 1e-11


def test_predefined_roots_output():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)

    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12,
                  method='adams', roots=roots, nroots=1)

    yout, info = integrate_predefined(f, None, [1], [0, 2], **kwargs)
    assert len(info['root_indices']) == 1
    roots_x, roots_y = info['roots_output']
    assert roots_x.shape == (1,)
    assert roots_y.shape == (1, 1)
    assert abs(roots_x[-1] - 1) < 1e-11
    assert abs(roots_y[-1, 0] - exp(1)) < 1e-11
