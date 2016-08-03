#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters.smarts

Module defines SMARTS filters.
"""

from rdkit import RDConfig
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

        >>> data = [
        ...         skchem.Mol.from_smiles('CC', name='ethane'),
        ...         skchem.Mol.from_smiles('c1ccccc1', name='benzene'),
        ...         skchem.Mol.from_smiles('c1ccccc1-c2c(C=O)ccnc2', name='big')
        ... ]

        >>> f = skchem.filters.SMARTSFilter({'benzene': 'c1ccccc1', 'pyridine': 'c1ccccn1', 'acetyl': 'C=O'}, agg='any')
        >>> f.transform(data, agg=False)
                acetyl benzene pyridine
        ethane   False   False    False
        benzene  False    True    False
        big       True    True     True

        >>> f.transform(data)
        ethane     False
        benzene     True
        big         True
        dtype: bool

        >>> f.filter(data)
        benzene                <Mol: c1ccccc1>
        big        <Mol: O=Cc1ccncc1-c1ccccc1>
        Name: structure, dtype: object

        >>> f.agg = all
        >>> f.filter(data)
        big    <Mol: O=Cc1ccncc1-c1ccccc1>
        Name: structure, dtype: object
    """

    def __init__(self, smarts, **kwargs):

        def read_smarts(s):
            if isinstance(s, str):
                return Mol.from_smarts(s, mergeHs=True)
            else:
                return s

        self.smarts = pd.Series(smarts).apply(read_smarts)
        super(SMARTSFilter, self).__init__(**kwargs)

    def _transform_mol(self, mol):

        return self.smarts.apply(lambda smarts: smarts in mol).values

    @property
    def columns(self):
        return self.smarts.index


class PAINSFilter(SMARTSFilter):

    """ Whether a molecule passes the Pan Assay INterference (PAINS) filters.

    These are supplied with RDKit, and were originally proposed by Baell et al.

    References:
        [The original paper](http://dx.doi.org/10.1021/jm901137j)

    Examples:

        Basic usage as a function on molecules:

        >>> import skchem
        >>> benzene = skchem.Mol.from_smiles('c1ccccc1', name='benzene')
        >>> pf = skchem.filters.PAINSFilter()
        >>> pf.transform(benzene)
        True
        >>> catechol = skchem.Mol.from_smiles('Oc1c(O)cccc1', name='catechol')
        >>> pf.transform(catechol)
        False

        >>> res = pf.transform(catechol, agg=False)
        >>> res[res]
        names
        catechol_A(92)    True
        Name: PAINSFilter, dtype: bool

        More useful in combination with pandas DataFrames:

        >>> data = [benzene, catechol]
        >>> pf.transform(data)
        benzene      True
        catechol    False
        dtype: bool

        >>> pf.filter(data)
        benzene    <Mol: c1ccccc1>
        Name: structure, dtype: object
    """

    def __init__(self):

        super(PAINSFilter, self).__init__(self._load_pains(), agg='not any')

    def _load_pains(cls):

        """ Load PAINS included in rdkit into a pandas dataframe and cache as class attribute. """

        if not hasattr(cls, '_pains'):
            path = os.path.join(RDConfig.RDDataDir, 'Pains', 'wehi_pains.csv')
            pains = pd.read_csv(path, names=['pains', 'names'])
            pains['names'] = pains.names.str.lstrip('<regId=').str.rstrip('>')
            pains = pains.set_index('names').pains.apply(Mol.from_smarts, mergeHs=True)
            cls._pains = pains
        return cls._pains
