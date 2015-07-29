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
    >>> import skchem
    >>> from skchem.descriptors import skchemize
    >>> from skchem.core import Mol
    >>> f = skchemize(Chem.RDKFingerprint)
    >>> m = Mol.from_smiles('c1ccccc1')
    >>> f(m)
    0       0
    1       0
    2       0
    3       0
    4       0
    5       0
    6       0
    7       0
    8       0
    9       0
    10      0
    11      0
    12      0
    13      0
    14      0
    15      0
    16      0
    17      0
    18      0
    19      0
    20      0
    21      0
    22      0
    23      0
    24      0
    25      0
    26      0
    27      0
    28      0
    29      0
           ..
    2018    0
    2019    0
    2020    0
    2021    0
    2022    0
    2023    0
    2024    0
    2025    0
    2026    0
    2027    0
    2028    0
    2029    0
    2030    0
    2031    0
    2032    0
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
    dtype: int64
    >>> from skchem.data import resource
    >>> df = skchem.read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
    >>> df.structure.apply(f)
          0     1     2     3     4     5     6     7     8     9     ...   2038  \\
    name                                                              ...          
    297      0     0     0     0     0     0     0     0     0     0  ...      0   
    6324     0     0     0     0     0     0     0     0     0     0  ...      0   
    6334     0     0     0     0     0     0     0     0     0     0  ...      0   
    <BLANKLINE>
          2039  2040  2041  2042  2043  2044  2045  2046  2047  
    name                                                        
    297      0     0     0     0     0     0     0     0     0  
    6324     0     0     0     0     0     0     0     0     0  
    6334     0     0     0     0     0     0     0     0     0  
    <BLANKLINE>
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

        if isinstance(obj, skc.core.Mol):
            return self._calculate_m(obj)

        elif isinstance(obj, pd.DataFrame) or isinstance(obj, pd.Series):
            return self._calculate_df(obj)

        else:
            raise NotImplementedError

    def _calculate_m(self, m):

        """ calculate the fingerprint for a molecule. """

        return self.func(m)

    def _calculate_df(self, df):

        """ calculate a fingerprint for a dataframe of molecules """

        return df.structure.apply(self.func)
