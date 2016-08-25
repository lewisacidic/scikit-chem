#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
## skchem.core.atom

Defining atoms in scikit-chem.
"""

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem.rdchem import GetPeriodicTable
from rdkit.Chem.AtomPairs.Utils import NumPiElectrons

from .base import ChemicalObject, PropertyView, ChemicalObjectView
from ..resource import PERIODIC_TABLE


RD_PT = GetPeriodicTable()


class Atom(Chem.rdchem.Atom, ChemicalObject):

    """ Object representing an Atom in scikit-chem. """

    @property
    def owner(self):

        """ skchem.Mol: the owning molecule.

        Warnings:
            This will seg fault if the atom is created manually.
        """
        from .mol import Mol
        return Mol.from_super(self.GetOwningMol())

    @property
    def bonds(self):

        """ tuple<skchem.Bonds>: the bonds to this `Atom`. """

        from .bond import Bond  # as bond imports atom, have to do it here
        return tuple(Bond.from_super(bond) for bond in self.GetBonds())

    def neighbours(self):

        """ tuple<Atom>: the neighbours of the atom. """

        return tuple(Atom.from_super(neigh) for neigh in self.GetNeighbors())

    @property
    def symbol(self):

        """ str: the element symbol of the atom. """

        return self.GetSymbol()

    @property
    def atomic_number(self):

        """ int: the atomic number of the atom. """

        return self.GetAtomicNum()

    @property
    def atomic_mass(self):

        """ float: the atomic mass of the atom in u. """

        return self.GetMass()

    @property
    def formal_charge(self):

        """ int: the formal charge. """

        return int(self.GetFormalCharge())

    @property
    def degree(self):

        """ int: the degree of the atom. """

        return self.GetDegree()

    @property
    def depleted_degree(self):

        """ int: the degree of the atom in the h depleted molecular graph. """

        return self.degree - self.n_instanced_hs

    @property
    def full_degree(self):

        """ int: the full degree of the atom in the h full molecular graph. """

        return self.degree + self.n_hs

    @property
    def valence_degree(self):

        """ int: the valence degree.

        $$ \delta_i^v = Z_i^v - h_i $$

        Where $ Z_i^v $ is the number of valence electrons and $ h_i $ is the
        number of hydrogens.
        """

        res = self.n_val_electrons - self.n_hs

        # this seems drastic...
        #if self.principal_quantum_number > 2:
        #    res /= self.atomic_number - self.n_val_electrons - 1

        return res

    @property
    def n_hs(self):

        """ int: the instanced, implicit and explicit number of hydrogens """

        return self.GetTotalNumHs() + self.n_instanced_hs

    @property
    def n_implicit_hs(self):

        """ int: the number of implicit hydrogens. """

        return self.GetNumExplicitHs()

    @property
    def n_explicit_hs(self):

        """ int: the number of explicit hydrogens. """

        return self.GetNumImplicitHs()

    @property
    def n_instanced_hs(self):

        """ int: The number of instanced hydrogens. """

        return sum(a.atomic_number == 1
                   for a in self.neighbours())

    @property
    def n_total_hs(self):

        """ int: the total number of hydrogens (according to rdkit). """

        return self.GetTotalNumHs()

    @property
    def n_val_electrons(self):

        """ int: the number of valence electrons. """

        return RD_PT.GetNOuterElecs(self.atomic_number)

    @property
    def n_pi_electrons(self):

        """ int: the number of pi electrons. """

        return NumPiElectrons(self)

    @property
    def n_lone_pairs(self):

        """ int: the number of lone pairs. """

        return int(0.5 * (self.n_val_electrons - self.full_degree -
                          self.formal_charge - self.n_pi_electrons))

    @property
    def explicit_valence(self):

        """ int: the explicit valence. """

        return self.GetExplicitValence()

    @property
    def implicit_valence(self):

        """ int: the implicit valence. """

        return self.GetImplicitValence()

    @property
    def valence(self):

        """ int: the valence. """

        return self.GetTotalValence()

    @property
    def hybridization_state(self):

        """ str: the hybridization state. """

        return self.GetHybridization().name

    @property
    def chiral_tag(self):

        """ int: the chiral tag. """

        return self.GetChiralTag()

    @property
    def cahn_ingold_prelog(self):

        """ The Cahn Ingold Prelog chirality indicator. """
        try:
            return self.props['_CIPCode']
        except KeyError:
            if self.owner is not None:
                Chem.FindMolChiralCenters(self.owner)
        return self.props.get('_CIPCode', None)

    @property
    def is_terminal(self):
        """ bool: whether the atom is terminal. """

        return self.depleted_degree == 1

    @property
    def is_aromatic(self):

        """ bool: whether the atom is aromatic. """

        return self.GetIsAromatic()

    @property
    def is_in_ring(self):

        """ bool: whether the atom is in a ring. """

        return any(b.is_in_ring for b in self.bonds)

    @property
    def van_der_waals_radius(self):

        """ float: the Van der Waals radius in angstroms. """

        return PERIODIC_TABLE.van_der_waals_radius[self.atomic_number]

    @property
    def van_der_waals_volume(self):

        """ float: the van der waals volume in angstroms^3.

        $\frac{4}{3} \pi r_v^3 $ """

        return PERIODIC_TABLE.van_der_waals_volume[self.atomic_number]

    _cov_dict = {
        6: {'SP': 0.60, 'SP2': 0.67, 'SP3': 0.77},
        7: {'SP': 0.55, 'SP2': 0.62, 'SP3': 0.74},
        8: {'SP2': 0.62, 'SP3': 0.74},
        9: {'SP3': 0.72},
        15: {'SP2': 1.00, 'SP3': 1.10},
        16: {'SP2': 0.97, 'SP3': 1.04},
        17: {'SP3': 0.99},
        35: {'SP3': 1.14},
        55: {'SP3': 1.33}
    }

    @property
    def covalent_radius(self):

        """ float: the covalent radius in angstroms. """

        if self.atomic_number in self._cov_dict.keys():
            return self._cov_dict[self.atomic_number][self.hybridization_state]
        else:
            return PERIODIC_TABLE.covalent_radius[self.atomic_number]

    @property
    def ionisation_energy(self):

        """ float: the first ionisation energy in eV. """

        return PERIODIC_TABLE.first_ionisation_energy[self.atomic_number]

    @property
    def electron_affinity(self):

        """ float: the first electron affinity in eV. """

        return PERIODIC_TABLE.electron_affinity[self.atomic_number]

    @property
    def principal_quantum_number(self):

        """ int: the principle quantum number. """

        return np.digitize(self.atomic_number,
                           np.array([1, 3, 11, 19, 37, 55, 87, 121]))

    @property
    def polarisability(self):

        """ float: the atomic polarisability in 10^{-20} m^3. """

        return PERIODIC_TABLE.atomic_polarisability[self.atomic_number]

    @property
    def pauling_electronegativity(self):

        """ float: the pauling electronegativity on Pauling scale. """

        return PERIODIC_TABLE.pauling_electronegativity[self.atomic_number]

    @property
    def sanderson_electronegativity(self):

        """ float: the sanderson electronegativity on Pauling scale. """

        return PERIODIC_TABLE.sanderson_electronegativity[self.atomic_number]

    @property
    def kier_hall_electronegativity(self):

        """ float: the hall-keir electronegativity. """

        if self.atomic_number == 1:
            return -0.2
        else:
            # py2 compat
            return float(self.valence_degree - self.depleted_degree) / \
                   (self.principal_quantum_number ** 2)

    @property
    def mcgowan_parameter(self):

        """ float: the mcgowan volume parameter"""

        return PERIODIC_TABLE.mcgowan_parameter[self.atomic_number]

    @property
    def kier_hall_alpha_contrib(self):

        """ float: the covalent radius in angstroms. """

        return self.covalent_radius / 0.77 - 1  # 0.77 is sp3 C

    @property
    def intrinsic_state(self):

        """ float: the intrinsic state of the atom. """

        # py2compat
        return (float(2 / self.principal_quantum_number) ** 2 *
                self.valence_degree + 1) / self.depleted_degree

    @property
    def hexcode(self):

        """ The hexcode to use as a color for the atom. """

        return PERIODIC_TABLE.hexcode[self.atomic_number]

    @property
    def props(self):

        """ PropertyView: rdkit properties of the atom. """

        if not hasattr(self, '_props'):
            self._props = PropertyView(self)
        return PropertyView(self)

    def __repr__(self):

        return '<{klass} element="{symbol}" at {address}>'.format(
            klass=self.__class__.__name__,
            symbol=self.symbol,
            address=hex(id(self))
            )

    def __str__(self):

        return self.symbol


class AtomView(ChemicalObjectView):

    def __getitem__(self, index):
        res = super(AtomView, self).__getitem__(index)
        if res is None:
            if np.abs(index) >= len(self):
                msg = 'Index {} out of range for molecule with' \
                    '{} atoms.'.format(index, len(self))
                raise IndexError(msg)

            # index is negative, so adding gives desired indexing from back
            if index < 0:
                index += len(self)

            return Atom.from_super(self.owner.GetAtomWithIdx(index))

        else:
            return res

    def __len__(self):
        return self.owner.GetNumAtoms()

    @property
    def symbol(self):

        """ np.array<str>: the symbols of the atoms in view """

        return np.array([atom.symbol for atom in self], dtype=(np.str_, 3))

    @property
    def atomic_number(self):

        """ np.array<int>: the atomic number of the atoms in view """

        return np.array([atom.atomic_number for atom in self])

    @property
    def atomic_mass(self):

        """ np.array<float>: the atomic mass of the atoms in view """

        return np.array([atom.atomic_mass for atom in self])

    @property
    def formal_charge(self):

        """ np.array<int>: the formal charge on the atoms in view """

        return np.array([atom.formal_charge for atom in self])

    @property
    def degree(self):

        """ np.array<int>: the degree of the atoms in view, according to
        rdkit. """

        return np.array([atom.degree for atom in self])

    @property
    def depleted_degree(self):

        """ np.array<int>: the degree of the atoms in the view in the
        h-depleted molecular graph.  """

        return np.array([atom.depleted_degree for atom in self])

    @property
    def full_degree(self):

        """ np.array<int>: the degree of the atoms in the view in the
        h-filled molecular graph. """

        return np.array([atom.full_degree for atom in self])

    @property
    def valence_degree(self):

        """ np.array<int>: the valence degree of the atoms in the view."""

        return np.array([atom.valence_degree for atom in self])

    @property
    def n_hs(self):

        """ np.array<int>: the number of hydrogens bonded to atoms in view. """

        return np.array([atom.n_hs for atom in self])

    @property
    def n_implicit_hs(self):

        """ np.array<int>: the number of implicit hydrogens bonded to atoms
        in view, according to rdkit.  """

        return np.array([atom.n_implicit_hs for atom in self])

    @property
    def n_explicit_hs(self):

        """ np.array<int>: the number of explicit  hydrogens bonded to atoms
        in view, according to rdkit.  """

        return np.array([atom.n_explicit_hs for atom in self])

    @property
    def n_instanced_hs(self):

        """ np.array<int>: the number of instanced hydrogens bonded to atoms
        in view.

        In this case, instanced means the number hs explicitly initialized as
        atoms. """

        return np.array([atom.n_instanced_hs for atom in self])

    @property
    def n_total_hs(self):

        """ np.array<int>: the number of total hydrogens bonded to atoms in
        view, according to rdkit.  """

        return np.array([atom.n_total_hs for atom in self])

    @property
    def n_val_electrons(self):

        """ np.array<int>: the number of valence electrons bonded to atoms
        in view.  """

        return np.array([atom.n_val_electrons for atom in self])

    @property
    def n_pi_electrons(self):

        """ np.array<int>: the number of pi electrons on atoms in view. """

        return np.array([atom.n_pi_electrons for atom in self])

    @property
    def n_lone_pairs(self):

        """ np.array<int>: the number of lone pairs on atoms in view. """

        return (0.5 * (self.n_val_electrons - self.full_degree -
                       self.formal_charge - self.n_pi_electrons)).astype(int)

    @property
    def explicit_valence(self):

        """ np.array<int>: the explicit valence of the atoms in view.. """

        return np.array([atom.explicit_valence for atom in self])

    @property
    def implicit_valence(self):

        """ np.array<int>: the explicit valence of the atoms in view. """

        return np.array([atom.implicit_valence for atom in self])

    @property
    def valence(self):

        """ np.array<int>: the valence of the atoms in view. """

        return np.array([atom.valence for atom in self])

    @property
    def hybridization_state(self):

        """ np.array<str>: the hybridization state of the atoms in view.

         One of 'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'UNSPECIFIED', 'OTHER'"""

        return np.array([atom.hybridization_state for atom in self])

    @property
    def chiral_tag(self):

        """ np.array<str>: the chiral tag of the atoms in view. """

        return np.array([atom.chiral_tag for atom in self])

    @property
    def cahn_ingold_prelog(self):

        """ np.array<str>: the CIP string representation of atoms in view. """

        return np.array([atom.cahn_ingold_prelog for atom in self],
                        (np.str_, 4))

    @property
    def is_terminal(self):

        """ np.array<bool>: whether the atoms in the view are terminal. """

        return self.depleted_degree == 1

    @property
    def is_aromatic(self):

        """ np.array<bool>: whether the atoms in the view are aromatic. """

        return np.array([atom.is_aromatic for atom in self])

    @property
    def is_in_ring(self):

        """ np.array<bool>: whether the atoms in the view are in a ring."""

        return np.array([atom.is_in_ring for atom in self])

    @property
    def van_der_waals_radius(self):

        """ np.array<float>: the Van der Waals radius of the atoms in the
        view. """

        return PERIODIC_TABLE.van_der_waals_radius[self.atomic_number].values

    @property
    def van_der_waals_volume(self):

        """ np.array<float>: the Van der Waals volume of the atoms in the
        view. """

        return PERIODIC_TABLE.van_der_waals_volume[self.atomic_number].values

    @property
    def covalent_radius(self):

        """ np.array<float>: the covalent radius of the atoms in the view. """

        return np.array([atom.covalent_radius for atom in self])

    @property
    def ionisation_energy(self):

        """ np.array<float>: the first ionisation energy of the atoms in the
        view. """

        return PERIODIC_TABLE.first_ionisation_energy[
            self.atomic_number].values

    @property
    def electron_affinity(self):

        """ np.array<float>: the electron affinity of the atoms in the
        view. """

        return PERIODIC_TABLE.electron_affinity[self.atomic_number].values

    @property
    def principal_quantum_number(self):

        """ np.array<float>: the principal quantum number of the atoms in the
        view. """

        return np.digitize(self.atomic_number,
                           np.array([1, 3, 11, 19, 37, 55, 87, 121]))

    @property
    def polarisability(self):

        """ np.array<float>: the atomic polarisability of the atoms in the
        view. """

        return PERIODIC_TABLE.atomic_polarisability[self.atomic_number].values

    @property
    def pauling_electronegativity(self):

        """ np.array<float>: the pauling electronegativity of the atoms in the
        view. """

        return PERIODIC_TABLE.pauling_electronegativity[
            self.atomic_number].values

    @property
    def sanderson_electronegativity(self):

        """ np.array<float>: the sanderson electronegativity of the atoms in
        the view. """

        return PERIODIC_TABLE.sanderson_electronegativity[
            self.atomic_number].values

    @property
    def kier_hall_electronegativity(self):

        """ np.array<float>: the hall kier electronegativity of the atoms in
        the view."""

        return np.array([atom.kier_hall_electronegativity for atom in self])

    @property
    def mcgowan_parameter(self):

        """ np.array<float>: the mcgowan parameter of the atoms in the
        iew. """

        return PERIODIC_TABLE.mcgowan_parameter[self.atomic_number].values

    @property
    def kier_hall_alpha_contrib(self):

        """ np.array<float>: the contribution to the kier hall alpha for each
        atom in the view. """

        return self.covalent_radius / 0.77 - 1

    @property
    def intrinsic_state(self):

        """ np.ndarray<float>: the intrinsic state of the atoms in the
        view. """

        #py2 compat
        return ((2 / self.principal_quantum_number).astype(float) ** 2 *
                self.valence_degree + 1) / self.depleted_degree

    @property
    def hexcode(self):

        """ The hexcode to use as a color for the atoms in the view. """

        return PERIODIC_TABLE.hexcode[self.atomic_number].values
    @property
    def index(self):

        """ pd.Index: an index for the atoms in the view. """

        return pd.RangeIndex(len(self), name='atom_idx')
