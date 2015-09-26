# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np


def _check_callable(f, j, x0, y0):
    ny = len(y0)
    _fout = np.empty(ny)
    _ret = f(x0, y0, _fout)
    if _ret is not None:
        raise ValueError("f() must return None")

    if j is None:
        return  # Not all methods require a jacobian

    _jmat_out = np.empty((ny, ny))
    _dfdx_out = np.empty(ny)
    _ret = j(x0, y0, _jmat_out, _dfdx_out)
    if _ret is not None:
        raise ValueError("j() must return None")


def _check_indexing(f, j, x0, y0):
    ny = len(y0)
    _fout_short = np.empty(ny-1)
    try:
        f(x0, y0, _fout_short)
    except (IndexError, ValueError):
        pass
    else:
        raise ValueError("All elements in fout not assigned in f()")

    if j is None:
        return  # Not all methods require a jacobian

    _dfdx_out = np.empty(ny)
    _jmat_out_short = np.empty((ny, ny-1))
    try:
        j(x0, y0, _jmat_out_short, _dfdx_out)
    except (IndexError, ValueError):
        pass
    else:
        raise ValueError("All elements in Jout not assigned in j()")

    _jmat_out = np.empty((ny, ny))
    _dfdx_out_short = np.empty(ny-1)
    try:
        j(x0, y0, _jmat_out, _dfdx_out_short)
    except (IndexError, ValueError):
        pass
    else:
        raise ValueError("All elements in dfdx_out not assigned in j()")
