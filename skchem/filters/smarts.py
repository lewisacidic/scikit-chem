#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters.smarts

Module defines SMARTS filters.
"""

from rdkit import RDConfig
import pandas as pd
import os
import pandas as pd

from .base import Filter
from ..core import Mol

class SMARTSFilter(Filter):

    """ Filter a molecule based on smarts.

    Args:
        smarts (pd.Series):
            A series of SMARTS to use in the filter.
        agg (function):
            Option specifying the mode of the filter.

            - None : No filtering takes place
            - any: If any of the substructures are in molecule return True.
            - all: If all of the substructures are in molecule.

    Examples:

        >>> import skchem
        >>> m1 = skchem.Mol.from_smiles('CC')
        >>> m2 = skchem.Mol.from_smiles('c1ccccc1')
        >>> m3 = skchem.Mol.from_smiles('c1ccccc1-c2c(C=O)ccnc2')
        >>> ms = pd.Series({'ethane': m1, 'benzene': m2, 'big': m3})
        >>> f = skchem.filters.SMARTSFilter({'benzene': 'c1ccccc1', 'pyridine': 'c1ccccn1', 'acetyl': 'C=O'})
        >>> f.transform(ms)
                acetyl benzene pyridine
        benzene  False    True    False
        big       True    True     True
        ethane   False   False    False

        >>> f.filter(ms, agg=any)
        benzene                <Mol: c1ccccc1>
        big        <Mol: O=Cc1ccncc1-c1ccccc1>
        dtype: object

        >>> f.filter(ms, agg=all)
        big    <Mol: O=Cc1ccncc1-c1ccccc1>
        dtype: object
    """

    def __init__(self, smarts, **kwargs):

        def read_smarts(s):
            if isinstance(s, str):
                return Mol.from_smarts(s, mergeHs=True)
            else:
                return s

        self.smarts = pd.Series(smarts).apply(read_smarts)

        self.index = self.smarts.index
        super(SMARTSFilter, self).__init__(self.func, **kwargs)

    def func(self, mol):

        return self.smarts.apply(lambda smarts: smarts in mol)



class PAINSFilter(SMARTSFilter):

    """ Whether a molecule passes the Pan Assay INterference (PAINS) filters.

    These are supplied with RDKit, and were originally proposed by Baell et al.

    References:
        [The original paper](http://dx.doi.org/10.1021/jm901137j)

    Examples:

        Basic usage as a function on molecules:

        >>> import skchem
        >>> m1 = skchem.Mol.from_smiles('c1ccccc1', name='benzene')
        >>> no_pains = PAINSFilter()
        >>> no_pains(m1)
        True
        >>> m2 = skchem.Mol.from_smiles('Oc1c(O)cccc1', name='catechol')
        >>> no_pains(m2)
        False

        More useful in combination with pandas DataFrames:

        >>> import gzip
        >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
        >>> data = skchem.read_sdf(sdf)
        >>> no_pains.transform(data).value_counts()
        True     3855
        False     482
        dtype: int64

        >>> len(no_pains.filter(data))
        3855
    """

    def __init__(self):

        super(PAINSFilter, self).__init__(self._load_pains(), agg=any, neg=True)

    def _load_pains(cls):

        """ Load PAINS included in rdkit into a pandas dataframe and cache as class attribute. """

        if not hasattr(cls, '_pains'):
            path = os.path.join(RDConfig.RDDataDir, 'Pains', 'wehi_pains.csv')
            pains = pd.read_csv(path, names=['pains', 'names'])
            pains['names'] = pains.names.str.lstrip('<regId=').str.rstrip('>')
            pains = pains.set_index('names').pains.apply(Mol.from_smarts, mergeHs=True)
            cls._pains = pains
        return cls._pains
