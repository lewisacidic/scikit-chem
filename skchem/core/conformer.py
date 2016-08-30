#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.conformer

Defining conformers in scikit-chem.
"""

import warnings

import rdkit.Chem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdDistGeom import EmbedMolecule, EmbedMultipleConfs

import numpy as np

from .base import ChemicalObject, ChemicalObjectView


class Conformer(rdkit.Chem.rdchem.Conformer, ChemicalObject):

    """ Class representing a Conformer in scikit-chem. """

    @property
    def owner(self):

        """ skchem.Mol: the owning molecule. """
        from .mol import Mol
        return Mol.from_super(self.GetOwningMol())

    @property
    def positions(self):

        """ np.ndarray: the atom positions in the conformer.

        Note:
            This is a copy of the data, not the data itself.  You cannot
            allocate to a slice of this.
        """

        # cant slice this array sadly.

        return np.array([self.GetAtomPosition(i) for i in range(len(self))])

    @positions.setter
    def positions(self, val):

        assert val.shape[0] == len(self), 'Positions must be given only for ' \
                                          'atoms in the molecule.'
        assert 1 < val.shape[1] <= 3, 'Positions must be 2 or 3 dimensional.'

        if val.shape[1] == 2:
            val = np.hstack((val, np.zeros((len(val), 1))))

        self.Set3D(bool((val[:, 2] != 0).any()))

        return [self.SetAtomPosition(i, v) for i, v in enumerate(val)]

    @property
    def centre_of_mass(self):

        """ np.array: the centre of mass of the comformer. """

        atomic_mass = self.owner.atoms.atomic_mass
        return atomic_mass.dot(self.positions) / atomic_mass.sum()

    @property
    def geometric_centre(self):

        """ np.array: the geometric centre of the conformer. """

        return self.positions.mean(axis=0)

    def centre_representation(self, centre_of_mass=True):

        """ Centre representation to the center of mass.

        Args:
            centre_of_mass (bool):
                Whether to use the masses of atoms to calculate the centre of
                mass, or just use the mean position coordinate.

        Returns:
            Conformer
        """

        if centre_of_mass:
            self.positions -= self.centre_of_mass
        else:
            self.positions -= self.geometric_centre

    def _inertia_tensor(self):

        """ Calculate the inertia tensor. """

        mass = self.owner.atoms.atomic_mass
        pos = self.positions

        return np.array([[((pos[:, (i % 3) - 1] ** 2 +
                            pos[:, (j % 3) - 2] ** 2 if (i == j)
                            else - pos[:, i] * pos[:, j]) * mass).sum()
                          for i in range(3)] for j in range(3)])

    def align_with_principal_axes(self):

        """ Align the reference frame with the principal axes of inertia. """

        eig_val, eig_vects = np.linalg.eigh(self._inertia_tensor())
        self.positions = self.positions.dot(eig_vects)

    def canonicalize(self):

        """ Center the reference frame at the centre of mass and """

        self.centre_representation(centre_of_mass=True)
        self.align_with_principal_axes()

    @property
    def id(self):

        """ The ID of the conformer. """

        return self.GetId()

    @id.setter
    def id(self, value):

        warnings.warn('Setting the conformer value directly '
                      'may cause issues.', UserWarning)
        self.SetId(int(value))

    @property
    def is_3d(self):

        """ bool: whether the conformer is three dimensional. """

        return self.Is3D()

    @is_3d.setter
    def is_3d(self, val):

        pos = self.positions
        pos[:, 2] = 0
        self.positions = pos

    def __len__(self):

        return self.GetNumAtoms()

    def __repr__(self):
        return '<{klass} id="{id}" at {address}>'.format(
            klass=self.__class__.__name__,
            id=self.GetId(),
            address=hex(id(self)))


class ConformerView(ChemicalObjectView):

    def __getitem__(self, index):

        res = super(ConformerView, self).__getitem__(index)
        if res is None:
            try:
                return Conformer.from_super(
                    self.owner.GetConformer(int(index)))
            except ValueError:
                raise IndexError(
                        'Conformer id {} not available (choose one '
                        'of {}).'.format(index, tuple(self.id)))
        else:
            return res

    def __setitem__(self, key, value):

        assert isinstance(value, Conformer), 'Only `Conformer`s can be added.'
        assert len(value) == len(self.owner.atoms), \
            '`Conformers` must have the same number of atoms as the `Mol`.'
        self.owner.RemoveConformer(int(key))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            value.id = key
        self.owner.AddConformer(value)

    def __delitem__(self, key):

        self.owner.RemoveConformer(key)

    def __iter__(self):
        return ConformerIterator(self)

    def __len__(self):

        return self.owner.GetNumConformers()

    def append(self, value):

        assert isinstance(value, rdkit.Chem.Conformer)
        self[max(self.id) + 1] = value

    def append_2d(self, **kwargs):

        """ Append a 2D conformer. """

        kwargs['clearConfs'] = False
        Compute2DCoords(self.owner, **kwargs)

    def append_3d(self, n_conformers=1, **kwargs):

        """ Append (a) 3D conformer(s), roughly embedded but not optimized.

        Args:
            n_conformers (int):
                The number of conformers to append.
            Further kwargs are passed to `EmbedMultipleConfs`.

        """

        kwargs.setdefault('numConfs', n_conformers)
        kwargs['clearConfs'] = False
        EmbedMultipleConfs(self.owner, **kwargs)

    @property
    def positions(self):

        return np.array([conformer.positions for conformer in self])

    @property
    def is_3d(self):

        return np.array([conformer.is_3d for conformer in self])

    @property
    def id(self):

        return np.array([conf.GetId() for conf in self.owner.GetConformers()])


class ConformerIterator(object):

    """  Iterator for chemical object views.  """

    def __init__(self, view):
        """ Create an iterator from a chemical object view. """
        self.view = view
        self._ids = view.id
        self._current = 0
        self._high = len(view)

    def __next__(self):
        if self._current >= self._high:
            raise StopIteration
        else:
            self._current += 1
            return self.view[self._ids[self._current - 1]]

    # py2 compat
    next = __next__
