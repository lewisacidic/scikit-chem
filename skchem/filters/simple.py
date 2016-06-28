#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""

# skchem.filters.simple

Simple filters for compounds.

"""

from .base import Filter
import pandas as pd


class ElementFilter(Filter):

    """ Filter by elements.

    Args:
        elements: A list of elements to filter with.  If an element not in
        the list is found in a molecule, return False, else return True.
    """
    def __init__(self, elements, **kwargs):

        self.elements = elements

        super(ElementFilter, self).__init__(self.func)

    def func(self, mol):

        return all(atom.element in self.elements for atom in mol.atoms)


class OrganicFilter(ElementFilter):

    # TODO: rewrite the docs

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
            >>> is_organic = skchem.filters.OrganicFilter()
            >>> is_organic(m1)
            True
            >>> m2 = skchem.Mol.from_smiles('[cH-]1cccc1.[cH-]1cccc1.[Fe+2]', \
                                            name='ferrocene')
            >>> is_organic(m2)
            False

            More useful in combination with pandas data frames:

            >>> import gzip
            >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
            >>> data = skchem.read_sdf(sdf)
            >>> is_organic.apply(data).value_counts()
            True     4253
            False      84
            Name: structure, dtype: int64

            >>> len(is_organic.filter(data))
            4253
            >>> len(is_organic.filter(data, neg=True))
            84
    """

    elements = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']

    def __init__(self):
        super(OrganicFilter, self).__init__(self.elements)


def n_atoms(mol, above=2, below=75, include_hydrogens=False):

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

    """

    assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)

    n_a = len(mol.atoms)
    if include_hydrogens:
        n_a += sum(atom.GetNumImplicitHs() for atom in mol.atoms)

    return above <= n_a < below

class AtomNumberFilter(Filter):

    """Filter for whether the number of atoms in a molecule falls in a defined interval.

    ``above <= n_atoms < below``

    Args:
        above (int):
            The lower threshold number of atoms (exclusive).
            Defaults to None.
        below (int):
            The higher threshold number of atoms (inclusive).
            Defaults to None.

    Args:
        >>> import skchem
        >>> import gzip
        >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
        >>> data = skchem.read_sdf(sdf)
        >>> f_natom = skchem.filters.AtomNumberFilter(above=3, below=60)
        >>> f_natom.apply(data).value_counts()
        True     4306
        False      31
        Name: structure, dtype: int64

        >>> len(f_natom.filter(data))
        4306
        >>> len(f_natom.filter(data, neg=True))
        31
    """

    def __init__(self, above=3, below=60, include_hydrogens=False, **kwargs):

        assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
        self.above = above
        self.below = below
        self.include_hydrogens = include_hydrogens

        super(AtomNumberFilter, self).__init__(n_atoms, above=self.above,
                                below=self.below,
                                include_hydrogens=self.include_hydrogens,
                                **kwargs)


def mass(mol, above=10, below=900):

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
    """

    return above <= mol.mass < below


class MassFilter(Filter):
    """ Filter whether a the molecular weight of a molecule is lower than a threshold.

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

    Examples:

        >>> import skchem
        >>> import gzip
        >>> sdf = gzip.open(skchem.data.resource('ames_mutagenicity.sdf.gz'))
        >>> data = skchem.read_sdf(sdf)
        >>> f_mass = skchem.filters.MassFilter(above=10, below=900)
        >>> f_mass.apply(data).value_counts()
        True     4312
        False      25
        Name: structure, dtype: int64

        >>> len(f_mass.filter(data))
        4312
        >>> len(f_mass.filter(data, neg=True))
        25
    """

    def __init__(self, above=3, below=900, **kwargs):

        assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
        self.above = above
        self.below = below

        super(MassFilter, self).__init__(mass, above=self.above,
                                below=self.below, **kwargs)
