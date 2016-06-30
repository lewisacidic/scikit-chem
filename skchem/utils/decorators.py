#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
# skchem.utils.decorators

Decorators for skchem functions.
"""

import pandas as pd
from functools import wraps

from .. import core

def takes_mol_series(func):
    @wraps(func)
    @takes_pandas
    def inner(inp, *args, **kwargs):
        if isinstance(inp, pd.DataFrame):
            inp = inp.structure
        return func(inp, *args, **kwargs)
    return inner

def method_takes_mol_series(func):
    @wraps(func)
    def inner(self, *args, **kwargs):
        return takes_mol_series(lambda *args, **kwargs: func(self, *args, **kwargs))(*args, **kwargs)
    return inner

def takes_pandas(func):
    @wraps(func)
    def inner(inp, *args, **kwargs):
        single = False
        if isinstance(inp, core.Mol):
            single = True
            inp = pd.Series({inp.name: inp})
        elif isinstance(inp, (list, tuple)):
            inp = pd.Series(inp, index=[mol.name for mol in inp])
        elif isinstance(inp, dict):
            inp = pd.Series(inp)
        elif isinstance(inp, (pd.Series, pd.DataFrame)):
            pass
        else:
            raise NotImplementedError('{} cannot take object of type: {}'.format(func.__name__, type(inp)))
        res = func(inp, *args, **kwargs)
        if single:
            return res.iloc[0]
        else:
            return res
    return inner

def method_takes_pandas(func):
    @wraps(func)
    def inner(self, *args, **kwargs):
        return takes_pandas(lambda *args, **kwargs: func(self, *args, **kwargs))(*args, **kwargs)
    return inner
