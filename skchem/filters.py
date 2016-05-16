#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.filters

Chemical filters are defined.

"""

import os

from rdkit import RDConfig
import pandas as pd

from .core import Mol


def _load_pains():

    """ Load PAINS included in rdkit into a pandas dataframe """

    path = os.path.join(RDConfig.RDDataDir, 'Pains', 'wehi_pains.csv')
    pains = pd.read_csv(path, names=['pains', 'names'])
    pains['names'] = pains.names.str.lstrip('<regId=').str.rstrip('>')
    return pains.set_index('names').pains.apply(Mol.from_smarts, mergeHs=True)

PAINS = _load_pains()
ORGANIC = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']

def is_organic(mol):

    """ Whether a molecule is organic.

    For the purpose of this function, an organic molecule is defined as having
    atoms with elements only in the set H, B, C, N, O, F, P, S, Cl, Br, I.

    Args:
        mol (skchem.Mol):
            The molecule to be tested.

    Returns:
        bool:
            Whether the molecule is organic.

    Examples:

            Basic usage as a function on molecules:

            >>> import skchem
            >>> m1 = skchem.Mol.from_smiles('c1ccccc1', name='benzene')
            >>> skchem.filters.is_organic(m1)
            True
            >>> m2 = skchem.Mol.from_smiles('[cH-]1cccc1.[cH-]1cccc1.[Fe+2]', \
                                            name='ferrocene')
            >>> skchem.filters.is_organic(m2)
            False

            More useful in combination with pandas data frames:

            >>> import gzip
            >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
            >>> data = skchem.read_sdf(sdf)
            >>> data.structure.apply(skchem.filters.is_organic).value_counts()
            True     4253
            False      84
            Name: structure, dtype: int64
    """

    return all(atom.element in ORGANIC for atom in mol.atoms)


def no_pains(mol):

    """ Whether a molecule passes the Pan Assay INterference (PAINS) filters.

    These are supplied with RDKit, and were originally proposed by Baell et al.

    Args:
        mol: (skchem.Mol):
            The molecule to be tested.

    Returns:
        bool:
            Whether the molecule passes all the pains filters.

    References:
        [The original paper](http://dx.doi.org/10.1021/jm901137j)

    Examples:

            Basic usage as a function on molecules:

            >>> import skchem
            >>> m1 = skchem.Mol.from_smiles('c1ccccc1', name='benzene')
            >>> skchem.filters.no_pains(m1)
            True
            >>> m2 = skchem.Mol.from_smiles('Oc1c(O)cccc1', name='catechol')
            >>> skchem.filters.no_pains(m2)
            False

            More useful in combination with pandas data frames:

            >>> import gzip
            >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
            >>> data = skchem.read_sdf(sdf)
            >>> data.structure.apply(skchem.filters.no_pains).value_counts()
            True     3855
            False     482
            Name: structure, dtype: int64
    """

    return all(PAINS.apply(lambda pains: pains not in mol))


def n_atoms(mol, above=None, below=None, include_hydrogens=False):

    """ Whether the number of atoms in a molecule falls in a defined interval.

    ``above <= n_atoms < below``

    Args:
        mol: (skchem.Mol):
            The molecule to be tested.
        above (int):
            The lower threshold number of atoms (exclusive).
            Defaults to None.
        below (int):
            The higher threshold number of atoms (inclusive).
            Defaults to None.

    Returns:
        bool:
            Whether the molecule has more atoms than the threshold.

    Examples:

        Basic usage as a function on molecules:

        >>> import skchem
        >>> m = skchem.Mol.from_smiles('c1ccccc1') # benzene has 6 atoms.

        Lower threshold:

        >>> skchem.filters.n_atoms(m, above=3)
        True
        >>> skchem.filters.n_atoms(m, above=8)
        False

        Higher threshold:

        >>> skchem.filters.n_atoms(m, below=8)
        True
        >>> skchem.filters.n_atoms(m, below=3)
        False

        Bounds work like Python slices - inclusive lower, exclusive upper:

        >>> skchem.filters.n_atoms(m, above=6)
        True
        >>> skchem.filters.n_atoms(m, below=6)
        False

        Both can be used at once:

        >>> skchem.filters.n_atoms(m, above=3, below=8)
        True

        Can include hydrogens:

        >>> skchem.filters.n_atoms(m, above=3, below=8, include_hydrogens=True)
        False
        >>> skchem.filters.n_atoms(m, above=9, below=14, include_hydrogens=True)
        True

        More useful in combination with pandas data frames:

        >>> import gzip
        >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
        >>> data = skchem.read_sdf(sdf)
        >>> data.structure.apply(skchem.filters.n_atoms, above=5, below=50).value_counts()
        True     4211
        False     126
        Name: structure, dtype: int64

    """
    if not above:
        above = 0
    if not below:
        below = 1000000 # arbitrarily large number

    n_a = len(mol.atoms)
    if include_hydrogens:
        n_a += sum(a.GetNumImplicitHs() for a in mol.atoms)

    assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
    return above <= n_a < below


def mass(mol, above=None, below=None):

    """ Whether a the molecular weight of a molecule is lower than a threshold.

    ``above <= mass < below``

    Args:
        mol: (skchem.Mol):
            The molecule to be tested.
        above (float):
            The lower threshold on the mass.
            Defaults to None.
        below (float):
            The higher threshold on the mass.
            Defaults to None.

    Returns:
        bool:
            Whether the mass of the molecule is lower than the threshold.

    Examples:
        Basic usage as a function on molecules:

        >>> import skchem
        >>> m = skchem.Mol.from_smiles('c1ccccc1') # benzene has M_r = 78.
        >>> skchem.filters.mass(m, above=70)
        True
        >>> skchem.filters.mass(m, above=80)
        False
        >>> skchem.filters.mass(m, below=80)
        True
        >>> skchem.filters.mass(m, below=70)
        False
        >>> skchem.filters.mass(m, above=70, below=80)
        True

        More useful in combination with pandas data frames:

        >>> import gzip
        >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
        >>> data = skchem.read_sdf(sdf)
        >>> data.structure.apply(skchem.filters.mass, below=900).value_counts()
        True     4312
        False      25
        Name: structure, dtype: int64
    """

    if not above:
        above = 0
    if not below:
        below = 1000000

    assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
    return above <= mol.mass < below
