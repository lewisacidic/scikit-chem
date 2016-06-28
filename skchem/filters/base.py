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

        Filter results can be found with `apply`:
        >>> ethane = skchem.Mol.from_smiles('CC', name='ethane')
        >>> is_named.apply(ethane)
        True

        >>> anonymous = skchem.Mol.from_smiles('c1ccccc1')
        >>> is_named.apply(anonymous)
        False

        The filter can also be used as a function:
        >>> ethane = skchem.Mol.from_smiles('CC', name='ethane')
        >>> is_named.apply(ethane)
        True

        Apply can take a series or dataframe:
        >>> mols = pd.Series({'anonymous': anonymous, 'ethane': ethane})
        >>> is_named.apply(mols)
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

    def __init__(self, func,  agg=None, neg=None, **kwargs):

        self.func = func
        if agg == None:
            agg = lambda x: x
        self.agg = agg
        self.neg = neg
        self.kwargs = kwargs

    @method_takes_mol_series
    def transform(self, mols, agg=None, neg=None):

        """ Apply the function and return the boolean values. """

        if agg is None:
            agg = self.agg

        if neg is None:
            neg = self.neg

        res = mols.apply(self.func, **self.kwargs)

        if isinstance(res, pd.DataFrame):
            res = res.apply(agg, axis=1)

        if neg:
            res = ~res

        return res

    @method_takes_pandas
    def filter(self, X, y=None, agg=None, neg=None):

        """ Apply the function and return filtered values.

        Args:
            X (pd.Series or pd.DataFrame):
                The compound dataframe. """

        if neg is None:
            neg = self.neg

        if agg is None:
            agg = self.agg

        res = self.transform(X, neg=neg, agg=agg)

        if y is None:
            return X[res]
        else:
            return X[res], y[res]

    def __call__(self, *args, **kwargs):

        return self.apply(*args, **kwargs)
