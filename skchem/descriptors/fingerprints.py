#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.descriptors.fingerprints

Fingerprinting classes and associated functions are defined.

"""

import pandas as pd
from rdkit.Chem import DataStructs
import numpy as np
import skchem as skc

def skchemize(func, columns=None, *args, **kwargs):
    """

    transform an RDKit fingerprinting function to work well with pandas

    >>> from rdkit import Chem
    >>> from skchem import *
    >>> f = skchemize(Chem.RDKFingerprint)
    >>> m = Mol.from_smiles('c1ccccc1')
    >>> f(m)
    0     0
    1     0
    2     0
    3     0
    4     0
    5     0
    6     0
    7     0
    8     0
    9     0
    10    0
    11    0
    12    0
    13    0
    14    0
    ...
    2033    0
    2034    0
    2035    0
    2036    0
    2037    0
    2038    0
    2039    0
    2040    0
    2041    0
    2042    0
    2043    0
    2044    0
    2045    0
    2046    0
    2047    0
    Length: 2048, dtype: int64

    >>> df = skchem.read_sdf('skchem/tests/test_resources/hydrocarbons.sdf')
    >>> df.structure.apply(f)
         0     1     2     3     4     5     6     7     8     9     ...   2038  \
    Name                                                              ...
    297      0     0     0     0     0     0     0     0     0     0  ...      0
    6324     0     0     0     0     0     0     0     0     0     0  ...      0
    6334     0     0     0     0     0     0     0     0     0     0  ...      0

          2039  2040  2041  2042  2043  2044  2045  2046  2047
    Name
    297      0     0     0     0     0     0     0     0     0
    6324     0     0     0     0     0     0     0     0     0
    6334     0     0     0     0     0     0     0     0     0

    [3 rows x 2048 columns]

    """

    def func_wrapper(m):

        """ Function that wraps an rdkit function allowing it to produce dataframes. """

        arr = np.array(0)
        DataStructs.ConvertToNumpyArray(func(m, *args, **kwargs), arr)

        return pd.Series(arr, index=columns)

    return func_wrapper

class Fingerprinter(object):

    """ Fingerprinter class. """

    def __init__(self, columns, func):
        self.columns = columns
        self.func = func

    @classmethod
    def from_rdkit_func(cls, func, columns=None, *args, **kwargs):

        """ Construct a Fingerprinter from an rdkit function. """

        return Fingerprinter(columns, func=skchemize(func, columns=columns, *args, **kwargs))

    def __call__(self, obj):
        return self.calculate(obj)

    def calculate(self, obj):

        """ calculate the fingerprint for the given object. """

        if obj.__class__ is skc.core.Mol:
            return self._calculate_m(obj)

        elif obj.__class__ in [pd.DataFrame, pd.Series]:
            return self._calculate_df(obj)

    def _calculate_m(self, m):

        """ calculate the fingerprint for a molecule. """

        return self.func(m)

    def _calculate_df(self, df):

        """ calculate a fingerprint for a dataframe of molecules """

        return df.structure.apply(self.func)
