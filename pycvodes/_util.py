# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

valid_arg_combs = [{}, {"lband", "uband"}, {"nnz"}]


def _check_jac_type(**kwargs):
    nonnull_opts = dict((k, v) for k, v in kwargs.items() if v is not None)
    if any(map(set(nonnull_opts).__eq__, valid_arg_combs)):
        pass
    else:
        raise ValueError("Couldn't determine jacobian type from given non-default options: {}".format(nonnull_opts))


def _get_jmat_out(ny, lband=None, uband=None, nnz=None):
    if lband is None and nnz is None:
        # jmat_out, dfdx_out
        return np.empty((ny, ny)), np.empty(ny)
    elif nnz is None:
        # jmat_out, dfdx_out
        return np.empty((1 + lband + uband, ny)), np.empty(ny)
    else:
        # data, colptrs, rowvals
        return np.empty(nnz), np.empty(ny + 1), np.empty(nnz)


def _check_callable(f, j, x0, y0, lband=None, uband=None, nnz=None):
    ny = len(y0)
    _fout = np.empty(ny)
    _ret = f(x0, y0, _fout)
    if _ret is not None:
        raise ValueError("f() must return None")

    if j is None:
        return  # Not all methods require a jacobian

    args = _get_jmat_out(ny, lband=lband, uband=uband,
                         nnz=nnz)
    _ret = j(x0, y0, *args)
    if _ret is not None:
        raise ValueError("j() must return None")


def _get_jmat_out_short(ny, lband=None, uband=None, nnz=None):
    if lband is None and nnz is None:
        # jmat_out, dfdx_out
        return np.empty((ny, ny - 1)), np.empty(ny)
    elif nnz is None:
        # jmat_out, dfdx_out
        return np.empty((1 + lband + uband, ny - 1)), np.empty(ny)
    else:
        # data, colptrs, rowvals
        return np.empty(nnz - 1), np.empty(ny), np.empty(nnz-1)


def _check_indexing(f, j, x0, y0, lband=None, uband=None, nnz=None):
    ny = len(y0)
    _fout_short = np.empty(ny - 1)
    try:
        f(x0, y0, _fout_short)
    except (IndexError, ValueError):
        pass
    else:
        raise ValueError("All elements in fout not assigned in f()")

    if j is None:
        return  # Not all methods require a jacobian

    args = _get_jmat_out_short(ny, lband=lband, uband=uband, nnz=nnz)
    try:
        j(x0, y0, *args)
    except (IndexError, ValueError):
        pass
    else:
        if nnz:
            raise ValueError("In one of (data, colptrs, rowvals), not all elements assigned in j()")
        else:
            raise ValueError("In either jmat_out or dfdx_out, not all elements assigned in j()")
