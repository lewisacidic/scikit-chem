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
    def __init__(self, check_meso=True, n_jobs=1, verbose=True):

        """ Initialize a `ChiralFilter` object.

        Args:
            check_meso (bool):
                Whether to include a check for meso compounds.

            n_jobs (int):
                The number of procesess to run the filter in.

            verbose (bool):
                Whether to output a progress bar.
        """

        self.check_meso = check_meso
        super(ChiralFilter, self).__init__(agg='any', n_jobs=n_jobs,
                                           verbose=verbose)

    @property
    def columns(self):
        return pd.Index(['is_chiral'])

    @staticmethod
    def is_meso(mol):
        """ Determines whether the molecule is meso.

        Meso compounds have chiral centres, but has a mirror plane allowing
        superposition.

        Examples:
            >>> import skchem

            >>> cf = skchem.filters.ChiralFilter()

            >>> meso = skchem.Mol.from_smiles('F[C@H](Br)[C@H](Br)F')

            >>> cf.is_meso(meso)
            True

            >>> non_meso = skchem.Mol.from_smiles('F[C@H](Br)[C@@H](Br)F')
            >>> cf.is_meso(non_meso)
            False
        """
        with Suppressor():
            info = InchiInfo(mol.to_inchi())
            return info.get_sp3_stereo()['main']['non-isotopic'][2]

    def _transform_mol(self, mol):

        is_meso = self.is_meso(mol) if self.check_meso else False

        return any(atom.GetChiralTag() for atom in mol.atoms) and not is_meso
