#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.utils.helpers

Module providing helper functions for scikit-chem
"""

from collections import Iterable
from functools import wraps

import numpy as np
import pandas as pd


def optional_second_method(func):
    @wraps(func)
    def inner(self, arg, second_arg=None, **kwargs):
        res = func(self, arg, **kwargs)
        if second_arg is not None:
            return res, second_arg
        else:
            return res
    return inner


def iterable_to_series(mols):
    assert isinstance(mols, Iterable), 'Object must be iterable.'
    assert not isinstance(mols, str), 'Object cannot be a string.'
    if isinstance(mols, dict):
        return pd.Series(mols, name='structure')
    else:
        return pd.Series(mols, index=[mol.name if mol.name else i for i, mol in enumerate(mols)], name='structure')


def nanarray(shape):
    """ Produce an array of NaN in provided shape.

    Args:
        shape (tuple):
            The shape of the nan array to produce.

    Returns:
        np.array

    """
    return np.repeat(np.nan, np.prod(shape)).reshape(*shape)


def squeeze(data, axis=None):

    """ Squeeze dimension for length 1 arrays.

    Args:
        data (pd.Series or pd.DataFrame or pd.Panel):
            The pandas object to squeeze.
        axis (int or tuple):
            The axes along which to squeeze.

    Returns:
        pd.Series or pd.DataFrame
    """

    if axis is None:
        axis = range(len(data.axes))
    elif isinstance(axis, int):
        axis = (axis,)
    return data.iloc[tuple([0 if len(a) == 1 and i in axis else slice(None) for i, a in enumerate(data.axes)])]


class Defaults(object):

    def __init__(self, defaults):

        self.defaults = defaults

    def get(self, val):
        if isinstance(val, str):
           return self.defaults.get(val)
        else:
            return val

