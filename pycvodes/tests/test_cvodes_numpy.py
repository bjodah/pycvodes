# -*- coding: utf-8 -*-
from math import exp, pi
import os
import numpy as np
import pytest
import itertools as it

from pycvodes import (
    integrate_adaptive, integrate_predefined, requires_jac, get_include, iterative_linsols
)


def requires_klu(_test):
    from .._config import env
    return pytest.mark.skipif(env.get('NO_KLU', '0') == '1',
                              reason="sparse jacobian tests require KLU to be enabled")(_test)


def high_precision(_test):
    from .._config import env
    return pytest.mark.skipif(env['SUNDIALS_PRECISION'] == "single",
                              reason="test atol and/or rtol is too high for single precision")(_test)


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


def _get_f_j(k, with_sparse=False):
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

    _colptrs = [0, 2, 4, 5]
    _rowvals = [0, 1, 1, 2, 2]

    def j_sparse(t, y, data, colptrs, rowvals):
        data[:] = -k0, k0, -k1, k1, -k2
        colptrs[:] = _colptrs
        rowvals[:] = _rowvals

    if with_sparse:
        return f, j, j_sparse, len(_rowvals)
    else:
        return f, j


def _gravity_f_j_jtimes(g):
    """ RHS, jacobian, and jacobian-vector product for
        particle falling under influence of gravity"""

    def f(t, y, fout):
        x, v = y
        fout[0] = v
        fout[1] = -g

    def jac(t, y, jmat_out, dfdx_out=None, fy=None):
        jmat_out[0, 0] = 0
        jmat_out[0, 1] = 1
        jmat_out[1, 0] = 0
        jmat_out[0, 1] = 0

    def jtimes(v, Jv, t, y, fy, user_data=None, tmp=None):
        Jv[0] = v[1]
        Jv[1] = 0

    return f, jac, jtimes


def gravity_analytic(g, y0, tout):
    """ Analytic solution to the above at times tout"""
    x0, v0 = y0
    t0 = tout[0]
    return np.column_stack((x0 - 0.5 * g * (tout - t0)**2, v0 - g * (tout - t0)))


methods = [('adams', 1.8, False),
           ('adams', 1.8, True),
           ('bdf', 10, False),
           ('bdf', 10, True)]


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
    from .._config import env

    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j = _get_f_j(k)
    kwargs = {'method': method, 'dx0': 1e-10}
    if use_jac:
        if banded:
            jac_callbacks = [bandify(j, 1, 0), None]
            kwargs['lband'] = 1
            kwargs['uband'] = 0
        else:
            jac_callbacks = [j, None]
    else:
        jac_callbacks = [None]

    for j in jac_callbacks:
        xout = np.linspace(0, 3, 31)
        if env['SUNDIALS_PRECISION'] == "single":
            atol = [1e-4, 3e-5, 2e-5]
            rtol = 1e-4
        else:
            atol = [1e-8, 3e-9, 2e-9]
            rtol = 1e-8
        # Run twice to catch possible side-effects:
        yout, nfo = integrate_predefined(f, j, y0, xout, atol, rtol, **kwargs)
        yout, nfo = integrate_predefined(f, j, y0, xout, atol, rtol, **kwargs)
        yref = decay_get_Cref(k, y0, xout)
        assert np.allclose(yout, yref,
                           rtol=forgiveness*rtol,
                           atol=forgiveness*max(atol))
        assert np.allclose(nfo['atol'], atol) and np.isclose(nfo['rtol'], rtol)
        assert nfo['nfev'] > 0
        if method in requires_jac and j is not None:
            assert nfo['njev'] > 0
            if os.name == 'posix':
                assert nfo['time_jac'] >= 0


@high_precision
def test_integrate_adaptive_tstop0():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, x0=0, xend=3, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)

    xout, yout, info = integrate_adaptive(f, j, y0, x0=0, xend=0, **kwargs)
    assert xout == [0]
    assert np.allclose(yout, y0)


@pytest.mark.parametrize("method,forgiveness,banded", methods)
def test_integrate_adaptive(method, forgiveness, banded):
    from .._config import env

    use_jac = method in requires_jac
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    if env['SUNDIALS_PRECISION'] == "single":
        atol, rtol = 1e-4, 1e-4
    else:
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
    assert np.allclose(info['atol'], [atol]) and np.isclose(info['rtol'], rtol)

    with pytest.raises(RuntimeError) as excinfo:
        kw = kwargs.copy()
        kw['atol'], kw['rtol'] = 1e-36, 1e-36
        integrate_adaptive(f, j, y0, **kw)
    assert 'acc' in str(excinfo.value).lower()

    with pytest.raises(RuntimeError) as excinfo:
        integrate_adaptive(f, j, y0, nsteps=7, **kwargs)
    assert 'maximum' in str(excinfo.value).lower()
    assert '7' in str(excinfo.value).lower()


@high_precision
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

    with pytest.raises(RuntimeError) as excinfo:
        integrate_predefined(f, None, [1], [0, 1, 2], nsteps=7, **kwargs)
    assert 'too_much_work' in str(excinfo.value).lower()
    assert '7' in str(excinfo.value).lower()


@high_precision
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


@high_precision
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


@high_precision
def test_roots_adaptive():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    assert len(info['root_indices']) == 1
    assert info['n_root_evals'] > 10
    assert np.min(np.abs(xout - 1)) < 1e-11


@high_precision
def test_quads_adaptive():
    k = 0.7

    def f(t, y, fout):
        fout[0] = -0.7*y[0]

    def quads(t, y, out):
        out[0] = t*y[0]
        out[1] = y[0]**2

    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  quads=quads, nquads=2)
    A, t0, duration = 42, 0, 4
    t, yout, info = integrate_adaptive(f, None, [A, 0, 0], t0, t0 + duration, **kwargs)
    assert np.allclose(yout[:, 0], 42*np.exp(-k*t))
    q0 = A/k**2 + (-A*k**2*t - A*k)*np.exp(-k*t)/k**3
    q1 = (1.0/2.0)*A**2/k - 1.0/2.0*A**2*np.exp(-2*k*t)/k
    assert np.allclose(info['quads'][:, 0], q0)
    assert np.allclose(info['quads'][:, 1], q1)


@high_precision
def test_quads_predefined():
    k = 0.7

    def f(t, y, fout):
        fout[0] = -0.7*y[0]

    def quads(t, y, out):
        out[0] = t*y[0]
        out[1] = y[0]**2

    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  quads=quads, nquads=2)
    A, t0, duration = 42, 0, 4
    t = np.linspace(t0, duration, 37)
    yout, info = integrate_predefined(f, None, [A, 0, 0], t, **kwargs)
    assert np.allclose(yout[:, 0], 42*np.exp(-k*t))
    q0 = A/k**2 + (-A*k**2*t - A*k)*np.exp(-k*t)/k**3
    q1 = (1.0/2.0)*A**2/k - 1.0/2.0*A**2*np.exp(-2*k*t)/k
    assert np.allclose(info['quads'][:, 0], q0)
    assert np.allclose(info['quads'][:, 1], q1)


@high_precision
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


@high_precision
def test_roots_predefined_autorestart():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1, nsteps=3, autorestart=4)
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


@high_precision
def test_return_on_root():
    def f(t, y, fout):
        fout[0] = y[0]

    def roots(t, y, out):
        out[0] = y[0] - exp(1)
    kwargs = dict(dx0=1e-12, atol=1e-12, rtol=1e-12, method='adams',
                  roots=roots, nroots=1, return_on_root=True,
                  return_on_error=True)
    xout, yout, info = integrate_adaptive(f, None, [1], 0, 2, **kwargs)
    assert len(info['root_indices']) == 1
    assert abs(xout[-1] - 1) < 1e-11
    assert abs(yout[-1, 0] - exp(1)) < 1e-11
    assert info['success']


@high_precision
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


def test_rhs_unrecoverable_error():
    def f(t, y, fout):
        fout[0] = y[0]
        return -1
    kwargs = dict(dx0=0.0, atol=1e-8, rtol=1e-8, method='adams')
    with pytest.raises(RuntimeError):
        yout, info = integrate_predefined(f, None, [1], [0, 1, 2], **kwargs)


@high_precision
def test_rhs_recoverable_error():
    global idx
    idx = -2

    def f(t, y, fout):
        global idx
        fout[0] = y[0]
        idx = idx + 1
        return 1 if 0 < idx < 3 else 0

    kwargs = dict(dx0=0.0, atol=1e-8, rtol=1e-8, method='adams')
    yout, info = integrate_predefined(f, None, [1], [0, 1, 2], **kwargs)


@high_precision
def test_adaptive_return_on_error():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(x0=0, xend=3, dx0=1e-10, atol=atol, rtol=rtol,
                  method='bdf')
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, nsteps=7, return_on_error=True, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=10*rtol,
                       atol=10*atol)
    assert xout.size == 8
    assert 1e-6 < xout[-1] < 1
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] == False  # noqa
    assert xout[-1] < kwargs['xend']  # obviously not strict


@high_precision
def test_adaptive_autorestart():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(x0=0, xend=3, dx0=1e-10, atol=atol, rtol=rtol,
                  method='BDF', nsteps=62, return_on_error=True,
                  autorestart=10, autonomous_exprs=True)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=10*rtol,
                       atol=10*atol)
    assert xout[-1] > 1e-6
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success']
    assert xout[-1] == kwargs['xend']


@high_precision
def test_predefined_autorestart():
    k = k0, k1, k2 = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    x0, xend = 0, 3
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol,
                  method='BDF', nsteps=7,
                  autorestart=10)
    f, j = _get_f_j(k)
    xout = np.linspace(x0, xend)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref,
                       rtol=10*rtol,
                       atol=10*atol)
    assert xout[-1] > 1e-6
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success']
    assert xout[-1] == xend


@high_precision
def test_predefined_return_on_error():
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0., 0.]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol,
                  method='bdf', return_on_error=True, nsteps=7)
    f, j = _get_f_j(k)
    xout = np.logspace(-3, 1)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert np.allclose(yout[:info['nreached'], :], yref[:info['nreached'], :],
                       rtol=10*rtol,
                       atol=10*atol)
    assert 10 < info['nreached'] < 40
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is False


@high_precision
def test_dx0cb():
    k = 1e23, 3.0, 4.0
    y0 = [.7, .0, .0]
    x0, xend = 0, 5
    kwargs = dict(atol=1e-8, rtol=1e-8, method='bdf', dx0cb=lambda x, y: y[0]*1e-30)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, x0, xend, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, atol=40*kwargs['atol'], rtol=40*kwargs['rtol'])
    assert info['n_steps'] < 1000
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is True
    assert xout[-1] == xend


@high_precision
def test_dx_max_cb():
    k = 1e23, 3.0, 4.0
    y0 = [.7, .0, .0]
    x0, xend = 0, 5
    kwargs = dict(atol=1e-8, rtol=1e-8, method='bdf', dx_max_cb=lambda x, y: 1e-3, nsteps=xend*1050)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, x0, xend, **kwargs)
    yref = decay_get_Cref(k, y0, xout)
    assert np.allclose(yout, yref, atol=40*kwargs['atol'], rtol=40*kwargs['rtol'])
    assert info['n_steps'] > 1000
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success'] is True
    assert xout[-1] == xend


@high_precision
def test_predefined_ew_ele():
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0., 0.]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol,
                  method='bdf', ew_ele=True)
    f, j = _get_f_j(k)
    xout = np.logspace(-3, 1)
    yout, info = integrate_predefined(f, j, y0, xout, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success']
    abs_ew_ele = np.abs(np.prod(info['ew_ele'], axis=1))
    assert np.all(abs_ew_ele < 1)


@high_precision
def test_adaptive_ew_ele():
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0., 0.]
    atol, rtol = 1e-8, 1e-8
    kwargs = dict(dx0=1e-10, atol=atol, rtol=rtol,
                  method='bdf', ew_ele=True)
    f, j = _get_f_j(k)
    xout, yout, info = integrate_adaptive(f, j, y0, 0, 10, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
    assert yout.shape[0] == xout.size
    assert info['nfev'] > 0
    assert info['njev'] > 0
    assert info['success']
    abs_ew_ele = np.abs(np.prod(info['ew_ele'], axis=1))
    assert np.all(abs_ew_ele < 1)


jtimes_params = list(it.product(iterative_linsols, [True, False]))


@high_precision
@pytest.mark.parametrize("linear_solver,with_jac", jtimes_params)
def test_jtimes_adaptive(linear_solver, with_jac):
    g = 9.81
    y0 = [1000.0, 0.0]
    atol, rtol = 1e-8, 1e-8
    f, jac, jtimes = _gravity_f_j_jtimes(g)
    if not with_jac:
        jac = None
    kwargs = dict(atol=atol, rtol=rtol,
                  method='bdf', linear_solver=linear_solver,
                  jtimes=jtimes)
    tout, yout, info = integrate_adaptive(f, jac, y0, 0, 10, **kwargs)
    yref = gravity_analytic(g, y0, tout)
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
    assert info['success']
    assert info['njvev'] > 0
    if not with_jac:
        assert info['njev'] == 0


@high_precision
@pytest.mark.parametrize("linear_solver,with_jac", jtimes_params)
def test_jtimes_predefined(linear_solver, with_jac):
    g = 9.81
    y0 = [1000.0, 0.0]
    atol, rtol = 1e-8, 1e-8
    f, jac, jtimes = _gravity_f_j_jtimes(g)
    if not with_jac:
        jac = None
    kwargs = dict(atol=atol, rtol=rtol,
                  method='bdf', linear_solver=linear_solver,
                  jtimes=jtimes)
    tout = np.linspace(0, 10, 100)
    yout, info = integrate_predefined(f, jac, y0, tout, **kwargs)
    yref = gravity_analytic(g, y0, tout)
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
    assert info['success']
    assert info['njvev'] > 0
    if not with_jac:
        assert info['njev'] == 0


@high_precision
@requires_klu
def test_sparse_equals_dense():
    try:
        from scipy.sparse import csc_matrix
    except ImportError:
        pytest.skip("Need scipy.sparse to compare sparse/dense jacobian output")
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    f, j, j_sparse, nnz = _get_f_j(k, with_sparse=True)

    n = len(y0)
    J_dense = np.empty((n, n))
    data = np.empty(nnz)
    colptrs = np.empty(n + 1)
    rowvals = np.empty(nnz)
    j(0, y0, J_dense)
    j_sparse(0, y0, data, colptrs, rowvals)

    J_sparse = csc_matrix((data, rowvals, colptrs), shape=(n, n)).toarray()
    assert np.allclose(J_sparse, J_dense)


@high_precision
@requires_klu
def test_sparse_jac_adaptive():
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    atol, rtol = 1e-8, 1e-8
    f, _, j_sparse, nnz = _get_f_j(k, with_sparse=True)
    kwargs = dict(atol=atol, rtol=rtol, method='bdf',
                  linear_solver='klu', nnz=nnz)
    xout, yout, info = integrate_adaptive(f, j_sparse, y0, 0, 10, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert info['success']
    assert info['njev'] > 0
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)


@high_precision
@requires_klu
def test_sparse_jac_predefined():
    k = 2.0, 3.0, 4.0
    y0 = [0.7, 0.3, 0.5]
    xout = np.linspace(0, 3, 31)
    atol, rtol = 1e-8, 1e-8
    f, _, j_sparse, nnz = _get_f_j(k, with_sparse=True)
    kwargs = dict(method='bdf', linear_solver='klu', nnz=nnz)
    yout, info = integrate_predefined(f, j_sparse, y0, xout, 1.0e-8, 1.0e-8, **kwargs)
    yref = decay_get_Cref(k, y0, xout - xout[0])
    assert info['success']
    assert info['njev'] > 0
    assert np.allclose(yout, yref, rtol=10*rtol, atol=10*atol)
