# -*- coding: utf-8 -*-
"""
Python binding for cvodes from the sundials library.
"""

from __future__ import division, absolute_import

from ._cvodes_numpy import adaptive, predefined, requires_jac, steppers
from ._util import _check_callable, _check_indexing
from ._release import __version__


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)


def integrate_adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol,
                       dx_min=.0, dx_max=.0, nsteps=500, nderiv=0,
                       roots=None, nroots=0, return_on_root=False,
                       check_callable=False, check_indexing=False, **kwargs):
    """
    Integrates a system of ordinary differential equations.

    Parameters
    ----------
    rhs: callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac: callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0: array_like
        initial values of the dependent variables
    x0: float
        initial value of the independent variable
    xend: float
        stopping value for the independent variable
    dx0: float
        initial step-size
    atol: float
        absolute tolerance
    rtol: float
        relative tolerance
    dx_min: float
        minimum step (default: 0.0)
    dx_max: float
        maximum step (default: 0.0)
    nsteps: int
        maximum number of steps (default: 500)
    nderiv: int
        number of derivatives (default: 0)
    roots: callback (default: None)
        with signature roots(x, yarr[:ny], out[:nroots]) -> None
    nroots: int (default: 0)
        number of root functions in roots
    return_on_root: bool (defaul: False)
        exit early (on first found root)
    check_callable: bool (default: False)
        perform signature sanity checks on ``rhs`` and ``jac``
    check_indexing: bool (default: False)
        perform item setting sanity checks on ``rhs`` and ``jac``.
    \*\*kwargs:
         'method': str
            One of: 'adams' or 'bdf' (default: 'bdf')
         'iter_type': str (default: 'default')
            One of: 'default', 'functional', 'newton'
         'linear_solver': str (default: 'default')
            One of: 'default', 'dense', 'banded', 'gmres',
            'gmres_classic', 'bicgstab', 'tfqmr'

    Returns
    -------
    (xout, yout, info):
        xout: 1-dimensional array of values for the independent variable
        yout: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0)
        info: dictionary with information about the integration

    """
    # Sanity checks to reduce risk of having a segfault:
    lband, uband = kwargs.get('lband', None), kwargs.get('uband', None)
    if check_callable:
        _check_callable(rhs, jac, x0, y0, lband, uband)

    if check_indexing:
        _check_indexing(rhs, jac, x0, y0, lband, uband)

    return adaptive(rhs, jac, y0, x0, xend, dx0, atol, rtol, dx_min, dx_max,
                    nsteps, nderiv, roots=roots, nroots=nroots,
                    return_on_root=return_on_root, **kwargs)


def integrate_predefined(rhs, jac, y0, xout, dx0, atol, rtol,
                         dx_min=.0, dx_max=.0, nsteps=500, nderiv=0,
                         roots=None, nroots=0,
                         check_callable=False, check_indexing=False, **kwargs):
    """
    Integrates a system of ordinary differential equations.

    Parameters
    ----------
    rhs: callable
        Function with signature f(t, y, fout) which modifies fout *inplace*.
    jac: callable
        Function with signature j(t, y, jmat_out, dfdx_out) which modifies
        jmat_out and dfdx_out *inplace*.
    y0: array_like
        initial values of the dependent variables
    xout: array_like
        values of the independent variable
    dx0: float
        initial step-size
    atol: float
        absolute tolerance
    rtol: float
        relative tolerance
    dx_min: float
        minimum step (default: 0.0)
    dx_max: float
        maximum step (default: 0.0)
    nsteps: int
        maximum number of steps (default: 500)
    nderiv: int
        number of derivatives (default: 0)
    roots: callback (default: None)
        With signature ``roots(x, yarr[:ny], out[:nroots]) -> None``,
        see info['root_indices'], note that xout is unaffected.
    nroots: int (default: 0)
        number of root functions in roots
    check_callable: bool (default: False)
        perform signature sanity checks on ``rhs`` and ``jac``
    check_indexing: bool (default: False)
        perform item setting sanity checks on ``rhs`` and ``jac``.
    \*\*kwargs:
         'method': str
            One of: 'adams' or 'bdf' (default: 'bdf')
         'iter_type': str (default: 'default')
            One of: 'default', 'functional', 'newton'
         'linear_solver': str (default: 'default')
            One of: 'default', 'dense', 'banded', 'gmres',
            'gmres_classic', 'bicgstab', 'tfqmr'

    Returns
    -------
    (result, info):
        result: 2-dimensional array of the dependent variables (axis 1) for
            values corresponding to xout (axis 0)
        info: dictionary with information about the integration

    """
    # Sanity checks to reduce risk of having a segfault:
    lband, uband = kwargs.get('lband', None), kwargs.get('uband', None)
    if check_callable:
        _check_callable(rhs, jac, xout[0], y0, lband, uband)

    if check_indexing:
        _check_indexing(rhs, jac, xout[0], y0, lband, uband)

    return predefined(rhs, jac, y0, xout, dx0, atol, rtol, dx_min, dx_max,
                      nsteps, nderiv, roots=roots, nroots=nroots, **kwargs)
