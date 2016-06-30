#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters

Chemical filters are defined.

"""

import pandas as pd

from ..utils import method_takes_mol_series, method_takes_pandas

def _identity(x):
    return x

def _identity_meth(_, x):
    return x

class Filter(object):

    """ Filter base class

    Args:
        func (function: Mol => bool):
            The function to use to filter the arguments.
        agg (function: iterable<bool> => bool):
            The aggregation to use in the filter, for example `any` or `all`.

    Examples:

        >>> import skchem

        Initialize the filter with a function:
        >>> is_named = skchem.filters.Filter(lambda m: m.name is not None)

        Filter results can be found with `transform`:
        >>> ethane = skchem.Mol.from_smiles('CC', name='ethane')
        >>> is_named.transform(ethane)
        True

        >>> anonymous = skchem.Mol.from_smiles('c1ccccc1')
        >>> is_named.transform(anonymous)
        False

        The filter can also be used as a function:
        >>> ethane = skchem.Mol.from_smiles('CC', name='ethane')
        >>> is_named.transform(ethane)
        True

        Apply can take a series or dataframe:
        >>> mols = pd.Series({'anonymous': anonymous, 'ethane': ethane})
        >>> is_named.transform(mols)
        anonymous    False
        ethane        True
        dtype: bool

        Using `filter` will drop out molecules that fail the test:
        >>> is_named.filter(mols)
        ethane    <Mol: CC>
        dtype: object

        Only failed are retained with the `neg` keyword argument:
        >>> is_named.filter(mols, neg=True)
        anonymous    <Mol: c1ccccc1>
        dtype: object
    """

    _DEFAULT_AGG = _identity_meth
    _DEFAULT_IS_NEG = False

    def __init__(self, func,  agg=None, neg=False, **kwargs):

        self.func = func
        self.agg = self._get_agg(agg) if agg is not None else self._DEFAULT_AGG
        self.neg = neg
        self.kwargs = kwargs

    @property
    def agg(self):
        return self._agg

    @agg.setter
    def agg(self, val):
        self._agg = self._get_agg(val)

    def _get_agg(self, val):
        if val is True:
            return self._DEFAULT_AGG
        elif val is False:
            return _identity
        elif val is None:
            return self.agg
        else:
            return val

    @property
    def neg(self):
        return self._get_neg(self._neg)

    @neg.setter
    def neg(self, val):

        # xor
        self._neg = self._get_neg(val)


    def _get_neg(self, val):
        if self._DEFAULT_IS_NEG:
            return not val
        else:
            return val

    @method_takes_mol_series
    def transform(self, mols, agg=None, neg=None):

        """ Apply the function and return the boolean values. """

        agg = self._get_agg(agg)

        if neg is None:
            neg = self._neg
        else:
            neg = self._get_neg(neg)

        res = self._transform(mols)

        if isinstance(res, pd.DataFrame):
            res = res.apply(agg, axis=1)

        if neg:
            res = ~res

        return res

    def _transform(self, ser):
        return ser.apply(self.func, **self.kwargs)

    @method_takes_pandas
    def filter(self, X, y=None, agg=None, neg=None):

        """ Apply the function and return filtered values.

        Args:
            X (pd.Series or pd.DataFrame):
                The compound dataframe. """

        if neg is None:
            neg = self.neg

        agg = self._get_agg(agg)

        res = self.transform(X, neg=neg, agg=agg)

        if y is None:
            return X[res]
        else:
            return X[res], y[res]

    def __call__(self, *args, **kwargs):

        return self.transform(*args, **kwargs)
