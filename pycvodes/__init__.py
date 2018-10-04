# -*- coding: utf-8 -*-
"""
Python binding for cvodes from the sundials library.
"""
from __future__ import division, absolute_import

import numpy as np

from ._cvodes import adaptive, predefined, requires_jac, steppers, fpes, sundials_version
from ._util import _check_callable, _check_indexing
from ._release import __version__
from ._config import env as config


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)


def integrate_adaptive(rhs, jac, y0, x0, xend, atol, rtol, dx0=.0,
                       dx_min=.0, dx_max=.0, nsteps=500, method=None, nderiv=0,
                       roots=None, nroots=0, return_on_root=False,
                       check_callable=False, check_indexing=False, **kwargs):
    """ Integrates a system of ordinary differential equations.

    Solves the initial value problem (IVP) defined by the user supplied
    arguments. The solver chooses at what values of the independent variable
    results should be reported.

    Parameters
    ----------
    rhs : callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac : callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0 : array_like
        Initial values of the dependent variables.
    x0 : float
        Initial value of the independent variable.
    xend : float
        Stopping value for the independent variable.
    dx0 : float
        Initial step-size.
    atol : float
        Absolute tolerance.
    rtol : float
        Relative tolerance.
    dx_min : float
        Minimum step (default: 0.0).
    dx_max : float
        Maximum step (default: 0.0).
    nsteps : int
        Maximum number of steps (default: 500).
    method : str
        One of: 'adams' or 'bdf' (default: 'bdf')
    nderiv : int
        Number of derivatives (default: 0).
    roots : callback
        With signature ``roots(x, yarr[:ny], out[:nroots]) -> None``.
    nroots : int
        Number of root functions in roots.
    return_on_root : bool
        Exit early (on first found root).
    check_callable : bool
        perform signature sanity checks on ``rhs`` and ``jac``
    check_indexing : bool
        Perform item setting sanity checks on ``rhs`` and ``jac``.
    \*\*kwargs:
        'lband': int
            Number of lower bands.
            Indexing: ``banded[row_i - col_i + uband, col_i]``.
        'uband': int
            Number of upper bands.
            Indexing: ``banded[row_i - col_i + uband, col_i]``.
        'iter_type': str (default: 'default')
            One of: 'default', 'functional', 'newton'
        'linear_solver': str (default: 'default')
            One of: 'default', 'dense', 'banded', 'gmres',
            'gmres_classic', 'bicgstab', 'tfqmr'
        'return_on_error': bool
            Returns on error without raising an excpetion (with ``'success'==False``).
        'autorestart': int
            Useful for autonomous systems where conditions change during integration.
            Will restart the integration with ``x==0``. Maximum number of steps is then
            given by ``autorestart * nsteps``.
        'record_rhs_xvals' : bool
            When True: will return x values for rhs calls in ``info['rhs_xvals']``.
        'record_jac_xvals' : bool
            When True will return x values for jac calls in ``info['jac_xvals']``.
        'record_order' : bool
            When True will return used time stepper order in ``info['orders']``.
        'record_fpe' : bool
            When True will return observed floating point errors in ``info['fpes']``. (see ``fpes``)
        'record_steps' : bool
            When True will return stepsizes taken in ``info['steps']``.
        'dx0cb': callable
            Callback for calculating dx0 (make sure to pass ``dx0==0.0``) to enable.
            Signature: ``f(x, y[:]) -> float``.
        'dx_max_cb': callable
            Callback for calculating dx_max.
            Signature: ``f(x, y[:]) -> float``.
        'autonomous_exprs' bool
            Whether expressions contain the independent variable. If not, autorestart
            is allowed to shift the independent variable to zero at restart).
        'with_jtimes': bool
            Whether to use analytic Jacobian for iterative solves. (default: False)
        'ew_ele': bool
            Whether to return error_weights, estimated_local_errors in info dict.
        'constraints': array
            Per component constraints 0.0: no constraint, 1.0: >=0, -1.0: <=0, 2.0: >0.0, -2.0: <0.0.

    Returns
    -------
    (xout, yout, info):
        xout: 1-dimensional array of values for the independent variable
        yout: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0).
        info: Dictionary with information about the integration.

    """
    # Sanity checks to reduce risk of having a segfault:
    lband, uband = kwargs.get('lband', None), kwargs.get('uband', None)
    if check_callable:
        _check_callable(rhs, jac, x0, y0, lband, uband)

    if check_indexing:
        _check_indexing(rhs, jac, x0, y0, lband, uband)

    return adaptive(rhs, jac, np.ascontiguousarray(y0, dtype=np.float64), x0, xend,
                    atol, rtol, method or ('adams' if jac is None else 'bdf'),
                    nsteps, dx0, dx_min, dx_max, nderiv=nderiv, roots=roots, nroots=nroots,
                    return_on_root=return_on_root, **kwargs)


def integrate_predefined(rhs, jac, y0, xout, atol, rtol, dx0=.0,
                         dx_min=.0, dx_max=.0, nsteps=500, method=None,
                         nderiv=0, roots=None, nroots=0, check_callable=False,
                         check_indexing=False, **kwargs):
    """ Integrates a system of ordinary differential equations.

    Solves the initial value problem (IVP) defined by the user supplied
    arguments. The user chooses at what values of the independent variable
    results should be reported.

    Parameters
    ----------
    rhs : callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac : callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0 : array_like
        Initial values of the dependent variables.
    xout : array_like
        Values of the independent variable.
    dx0 : float
        Initial step-size.
    atol : float
        Absolute tolerance.
    rtol : float
        Relative tolerance.
    dx_min : float
        Minimum step (default: 0.0).
    dx_max : float
        Maximum step (default: 0.0).
    nsteps : int
        Maximum number of steps (default: 500).
    method : str
        One of: 'adams' or 'bdf' (default: 'bdf').
    nderiv : int
        Number of derivatives (default: 0).
    roots : callback (default: None)
        With signature ``roots(x, yarr[:ny], out[:nroots]) -> None``,
        see info['root_indices'], note that xout is unaffected.
    nroots : int (default: 0)
        Number of root functions in roots.
    check_callable : bool (default: False)
        Perform signature sanity checks on ``rhs`` and ``jac``.
    check_indexing : bool (default: False)
        Perform item setting sanity checks on ``rhs`` and ``jac``.
    \*\*kwargs:
        'lband': int
            Number of lower bands.
            Indexing: ``banded[row_i - col_i + uband, col_i]``.
        'uband': int
            Number of upper bands.
            Indexing: ``banded[row_i - col_i + uband, col_i]``.
        'iter_type': str (default: 'default')
            One of: 'default', 'functional', 'newton'.
        'linear_solver': str (default: 'default')
            One of: 'default', 'dense', 'banded', 'gmres',
            'gmres_classic', 'bicgstab', 'tfqmr'.
        'return_on_error': bool
            Returns on error without raising an excpetion (with ``'success'==False``).
        'autorestart': int
            Useful for autonomous systems where conditions change during integration.
            Will restart the integration with ``x==0``. Maximum number of steps is then
            given by ``2**autorestart * nsteps``.
        'record_rhs_xvals' : bool
            When True: will return x values for rhs calls in ``info['rhs_xvals']``.
        'record_jac_xvals' : bool
            When True will return x values for jac calls in ``info['jac_xvals']``.
        'record_order' : bool
            When True will return used time stepper order in ``info['orders']``.
        'record_fpe' : bool
            When True will return observed floating point errors in ``info['fpes']``. (see ``fpes``)
        'record_steps' : bool
            When True will return stepsizes taken in ``info['steps']``.
        'dx0cb': callable
            Callback for calculating dx0 (make sure to pass ``dx0==0.0``) to enable.
            Signature: ``f(x, y[:]) -> float``.
        'dx_max_cb': callable
            Callback for calculating dx_max.
            Signature: ``f(x, y[:]) -> float``.
        'autonomous_exprs' bool
            Whether expressions contain the independent variable. If not, autorestart
            is allowed to shift the independent variable to zero at restart).
        'with_jtimes': bool
            Whether to use analytic Jacobian for iterative solves. (default: False)
        'ew_ele': bool
            Whether to return error_weights, estimated_local_errors in info dict.
        'constraints': array
            Per component constraints 0.0: no constraint, 1.0: >=0, -1.0: <=0, 2.0: >0.0, -2.0: <0.0.

    Returns
    -------
    (yout, info):
        yout: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0)
        info: Dictionary with information about the integration.

    """
    # Sanity checks to reduce risk of having a segfault:
    lband, uband = kwargs.get('lband', None), kwargs.get('uband', None)
    if check_callable:
        _check_callable(rhs, jac, xout[0], y0, lband, uband)

    if check_indexing:
        _check_indexing(rhs, jac, xout[0], y0, lband, uband)

    return predefined(
        rhs, jac,
        np.ascontiguousarray(y0, dtype=np.float64),
        np.ascontiguousarray(xout, dtype=np.float64),
        atol, rtol, method or ('adams' if jac is None else 'bdf'),
        nsteps, dx0, dx_min, dx_max, nderiv=nderiv, roots=roots,
        nroots=nroots, **kwargs)
