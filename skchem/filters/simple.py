#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""

# skchem.filters.simple

Simple filters for compounds.

"""

from collections import Counter

import numpy as np
import pandas as pd

from ..resource import ORGANIC, PERIODIC_TABLE
from .base import Filter


class ElementFilter(Filter):

    """ Filter by elements.

        Args:
            elements (list[str]):
                A list of elements to filter with.  If an element not in the list is
                found in a molecule, return False, else return True.

            as_bits (bool):
                Whether to return integer counts or booleans for atoms if mode is `count`.

        Examples:

            Basic usage on molecules:

            >>> import skchem
            >>> has_halogen = skchem.filters.ElementFilter(['F', 'Cl', 'Br', 'I'], agg='any')

            Molecules with one of the atoms transform to `True`.

            >>> m1 = skchem.Mol.from_smiles('ClC(Cl)Cl', name='chloroform')
            >>> has_halogen.transform(m1)
            True

            Molecules with none of the atoms transform to `False`.

            >>> m2 = skchem.Mol.from_smiles('CC', name='ethane')
            >>> has_halogen.transform(m2)
            False

            Can see the atom breakdown by passing `agg` == `False`:
            >>> has_halogen.transform(m1, agg=False)
            has_element
            F     0
            Cl    3
            Br    0
            I     0
            Name: ElementFilter, dtype: int64

            Can transform series.

            >>> ms = [m1, m2]
            >>> has_halogen.transform(ms)
            chloroform     True
            ethane        False
            dtype: bool

            >>> has_halogen.transform(ms, agg=False)
            has_element  F  Cl  Br  I
            chloroform   0   3   0  0
            ethane       0   0   0  0

            Can also filter series:

            >>> has_halogen.filter(ms)
            chloroform    <Mol: ClC(Cl)Cl>
            Name: structure, dtype: object

            >>> has_halogen.filter(ms, neg=True)
            ethane    <Mol: CC>
            Name: structure, dtype: object

        """
    def __init__(self, elements=None, as_bits=False, **kwargs):
        self.elements = elements
        self.as_bits = as_bits
        super(ElementFilter, self).__init__(**kwargs)

    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, val):
        if val is None:
            self._elements = PERIODIC_TABLE.symbol.tolist()
        else:
            self._elements = val

    @property
    def columns(self):
        return pd.Index(self.elements, name='has_element')

    def _transform_mol(self, mol):

        counter = Counter(atom.element for atom in mol.atoms)
        res = pd.Series(counter)

        res = res[self.elements].fillna(0).astype(int)

        if self.as_bits:
            res = (res > 0).astype(np.uint8)

        return res


class OrganicFilter(ElementFilter):
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
                >>> of = skchem.filters.OrganicFilter()
                >>> benzene = skchem.Mol.from_smiles('c1ccccc1', name='benzene')

                >>> of.transform(benzene)
                True

                >>> ferrocene = skchem.Mol.from_smiles('[cH-]1cccc1.[cH-]1cccc1.[Fe+2]',
                ...                                    name='ferrocene')
                >>> of.transform(ferrocene)
                False

                More useful on collections:

                >>> sa = skchem.Mol.from_smiles('CC(=O)[O-].[Na+]', name='sodium acetate')
                >>> norbornane = skchem.Mol.from_smiles('C12CCC(C2)CC1', name='norbornane')

                >>> data = [benzene, ferrocene, norbornane, sa]
                >>> of.transform(data)
                benzene            True
                ferrocene         False
                norbornane         True
                sodium acetate    False
                dtype: bool

                >>> of.filter(data)
                benzene          <Mol: c1ccccc1>
                norbornane    <Mol: C1CC2CCC1C2>
                Name: structure, dtype: object

                >>> of.filter(data, neg=True)
                ferrocene         <Mol: [Fe+2].c1cc[cH-]c1.c1cc[cH-]c1>
                sodium acetate                  <Mol: CC(=O)[O-].[Na+]>
                Name: structure, dtype: object
        """

    def __init__(self):
        super(OrganicFilter, self).__init__(elements=None, agg='not any')
        self.elements = [element for element in self.elements if element not in ORGANIC]


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
        n_a += sum(atom.GetNumImplicitHs() + atom.GetNumExplicitHs() for atom in mol.atoms)

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

    Examples:
        >>> import skchem

        >>> data = [
        ...         skchem.Mol.from_smiles('CC', name='ethane'),
        ...         skchem.Mol.from_smiles('CCCC', name='butane'),
        ...         skchem.Mol.from_smiles('NC(C)C(=O)O', name='alanine'),
        ...         skchem.Mol.from_smiles('C12C=CC(C=C2)C=C1', name='barrelene')
        ... ]

        >>> af = skchem.filters.AtomNumberFilter(above=3, below=7)

        >>> af.transform(data)
        ethane       False
        butane        True
        alanine       True
        barrelene    False
        Name: num_atoms_in_range, dtype: bool

        >>> af.filter(data)
        butane            <Mol: CCCC>
        alanine    <Mol: CC(N)C(=O)O>
        Name: structure, dtype: object

        >>> af = skchem.filters.AtomNumberFilter(above=5, below=15, include_hydrogens=True)

        >>> af.transform(data)
        ethane        True
        butane        True
        alanine       True
        barrelene    False
        Name: num_atoms_in_range, dtype: bool
    """

    def __init__(self, above=3, below=60, include_hydrogens=False, **kwargs):

        assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
        self.above = above
        self.below = below
        self.include_hydrogens = include_hydrogens

        super(AtomNumberFilter, self).__init__(**kwargs)

    def _transform_mol(self, mol):
        return n_atoms(mol, above=self.above, below=self.below, include_hydrogens=self.include_hydrogens)

    @property
    def columns(self):
        return pd.Index(['num_atoms_in_range'])

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

        >>> data = [
        ...         skchem.Mol.from_smiles('CC', name='ethane'),
        ...         skchem.Mol.from_smiles('CCCC', name='butane'),
        ...         skchem.Mol.from_smiles('NC(C)C(=O)O', name='alanine'),
        ...         skchem.Mol.from_smiles('C12C=CC(C=C2)C=C1', name='barrelene')
        ... ]

        >>> mf = skchem.filters.MassFilter(above=31, below=100)

        >>> mf.transform(data)
        ethane       False
        butane        True
        alanine       True
        barrelene    False
        Name: mass_in_range, dtype: bool

        >>> mf.filter(data)
        butane            <Mol: CCCC>
        alanine    <Mol: CC(N)C(=O)O>
        Name: structure, dtype: object

    """

    def __init__(self, above=3, below=900, **kwargs):

        assert above < below, 'Interval {} < a < {} undefined.'.format(above, below)
        self.above = above
        self.below = below

        super(MassFilter, self).__init__( **kwargs)

    def _transform_mol(self, mol):
        return mass(mol, above=self.above, below=self.below)

    @property
    def columns(self):
        return pd.Index(['mass_in_range'])
