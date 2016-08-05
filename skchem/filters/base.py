#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters

Chemical filters are defined.

"""

import pandas as pd

from ..base import BaseTransformer, Transformer
from .. import core
from ..utils import iterable_to_series, Defaults, optional_second_method


AGGS = Defaults(defaults={
    'none': lambda x: x,
    'any': any,
    'all': all,
    'not all': lambda x: not all(x),
    'not any': lambda x: not any(x)
})


class BaseFilter(BaseTransformer):

    def __init__(self, agg='any', **kwargs):
        super(BaseFilter, self).__init__(**kwargs)
        self.agg = agg

    @property
    def agg(self):
        """ callable: The aggregate function to use.  String aliases
        for `'any'`, `'not any'`, 'all', `'not all'` are available."""
        return self._agg

    @agg.setter
    def agg(self, val):
        self._agg = AGGS.get(val)

    @property
    def columns(self):

        """ pd.Index: The column index to use. """

        return pd.Index([self.__class__.__name__])

    def _mask(self, mols=None, res=None, neg=False):

        """ Generate a mask from molecules, or from their result after transform.

        Args:
            mols (pd.Series<skchem.Mol>):
                The molecules to use to generate the mask.
            res (pd.Series):
                The result of a transform. Overrides mols.
            neg (bool):
                Whether the mask should be inversed.

        Returns:
            pd.Series<bool>
        """

        res = self.transform(mols, agg=False) if res is None else res
        res = (res != False) & pd.notnull(res)
        if isinstance(res, pd.Series) and isinstance(mols, core.Mol):
            res = self.agg(res)
        if isinstance(res, pd.DataFrame):
            res = res.apply(self.agg, axis=1)
        return res == False if neg else res

    @optional_second_method
    def transform(self, mols, agg=True, **kwargs):

        # transform takes additional optional kwarg `agg`, that specifies to transform to the aggregated value or
        # return the full series.

        if agg:
            return self._mask(mols)
        else:
            return super(BaseFilter, self).transform(mols, **kwargs)

    def filter(self, mols, y=None, neg=False):

        mask = self._mask(mols=mols, neg=neg)

        if isinstance(mols, core.Mol):
            return mols if mask else None

        elif not isinstance(mols, pd.Series):
            mols = iterable_to_series(mols)

        if y is None:
            return mols[mask]
        else:
            return mols[mask], y[mask]


class Filter(BaseFilter, Transformer):
    """ Filter base class.

     Args:
         func (function: Mol => bool):
             The function to use to filter the arguments.
         agg (str or function: iterable<bool> => bool):
             The aggregation to use in the filter. Can be 'any', 'all', 'not any', 'not all' or a callable,
             for example `any` or `all`.

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

         Can take a series or dataframe:
         >>> mols = pd.Series({'anonymous': anonymous, 'ethane': ethane})
         >>> is_named.transform(mols)
         anonymous    False
         ethane        True
         Name: Filter, dtype: bool

         Using `filter` will drop out molecules that fail the test:
         >>> is_named.filter(mols)
         ethane    <Mol: CC>
         dtype: object

         Only failed are retained with the `neg` keyword argument:
         >>> is_named.filter(mols, neg=True)
         anonymous    <Mol: c1ccccc1>
         dtype: object
     """
    def __init__(self, func=None, **kwargs):
        super(Filter, self).__init__(**kwargs)
        if func is not None:
            self._transform_mol = func

    def _transform_mol(self, mol):
        raise NotImplemented


class TransformFilter(BaseFilter):

    """ Transform Filter object.

     Implements `transform_filter`, which allows a transform, then a
     filter step returning the transformed values that are not `False`, `None` or `np.nan`.

     """

    def transform_filter(self, mols, y=None, neg=False):

        res = self.transform(mols)
        mask = self._mask(res=res, neg=neg)

        if isinstance(mols, core.Mol):
            return res if mask else None

        if y is None:
            return res[mask]
        else:
            return res[mask], y[mask]