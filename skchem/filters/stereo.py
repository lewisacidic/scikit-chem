#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters.stereo

Stereo filters for scikit-chem.
"""

from rdkit.Chem.MolKey.InchiInfo import InchiInfo
import pandas as pd

from .base import Filter
from ..utils import Suppressor

class ChiralFilter(Filter):

    """ Filter chiral compounds.

    Examples:
        >>> import skchem
        >>> cf = skchem.filters.ChiralFilter()
        >>> ms = [
        ...    skchem.Mol.from_smiles('F[C@@H](F)[C@H](F)F', name='achiral'),
        ...    skchem.Mol.from_smiles('F[C@@H](Br)[C@H](Br)F', name='chiral'),
        ...    skchem.Mol.from_smiles('F[C@H](Br)[C@H](Br)F', name='meso'),
        ...    skchem.Mol.from_smiles('FC(Br)C(Br)F', name='racemic')
        ... ]
        >>> cf.transform(ms)
        achiral    False
        chiral      True
        meso       False
        racemic    False
        Name: is_chiral, dtype: bool

    """
    def __init__(self, check_meso=True, **kwargs):
        self.check_meso = True
        super(ChiralFilter, self).__init__(**kwargs)

    @property
    def columns(self):
        return pd.Index(['is_chiral'])

    def _is_meso(self, mol):
        """ Determines whether the molecule is meso (i.e. has chiral centres, but has a mirror plane allowing
        superposition).

        Examples:
            >>> import skchem

            >>> cf = skchem.filters.ChiralFilter()

            >>> meso = skchem.Mol.from_smiles('F[C@H](Br)[C@H](Br)F')

            >>> cf._is_meso(meso)
            True

            >>> non_meso = skchem.Mol.from_smiles('F[C@H](Br)[C@@H](Br)F')
            >>> cf._is_meso(non_meso)
            False
        """
        with Suppressor():
            return InchiInfo(mol.to_inchi()).get_sp3_stereo()['main']['non-isotopic'][2]

    def _transform_mol(self, mol):

        is_meso = self._is_meso(mol) if self.check_meso else False

        return any(atom.GetChiralTag() for atom in mol.atoms) and not is_meso
