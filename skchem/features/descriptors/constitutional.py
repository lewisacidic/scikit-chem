#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.features.descriptors.constitutional

Constitutional features for scikit-chem.
"""

from collections import OrderedDict
from functools import partial

from rdkit.Chem import rdMolDescriptors, rdmolops, Descriptors

from .caching import cache, requires_h_filled, requires_h_depleted
from .fundamentals import bond_order_adjacency_matrix


def molecular_weight(mol):

    """ The molecular weight.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return rdMolDescriptors.CalcExactMolWt(mol)


@requires_h_filled
def average_molecular_weight(mol):

    """ Average molecular weight of atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return molecular_weight(mol) / len(mol.atoms)


def sum_van_der_waals_volume(mol):

    """ The sum of the Van der Waals volume.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.van_der_waals_volume.sum()


def sum_electronegativity(mol):

    """ The sum of the Sanderson electronegativities.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.sanderson_electronegativity.sum()


def sum_ionisation_energy(mol):

    """ The sum of the first ionisation energies.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.ionisation_energy.sum()


def sum_polarisability(mol):

    """ The sum of the polarisabilities.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.polarisability.sum()


def mean_van_der_waals_volume(mol):

    """ The mean of the Van der Waals volume.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.van_der_waals_volume.sum()


def mean_electronegativity(mol):

    """ The mean of the Sanderson electronegativity.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.sanderson_electronegativity.mean()


def mean_ionisation_energy(mol):

    """ The mean of the first ionisation energies.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.ionisation_energy.mean()


def mean_polarisability(mol):

    """ The mean of the polarizabilities.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return mol.atoms.polarisability.mean()


@requires_h_depleted
def graph_density(mol):

    """ The graph density of the h-depleted graph.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return 2 * len(mol.bonds) / (len(mol.atoms) * (len(mol.atoms) - 1))


def n_atoms(mol):

    """ The number of atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return len(mol.atoms)


@requires_h_filled
def n_hyd(mol):

    """ The number of hydrogen atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return sum(mol.atoms.atomic_number == 1)


@requires_h_depleted
def n_atom(mol, symbol='C'):

    """ The number of atoms of symbol *s*.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return sum(mol.atoms.symbol == symbol)


@requires_h_filled
def fract_atom(mol, symbol='C'):

    """ The fraction of atoms of symbol *s*.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return sum(mol.atoms.symbol == symbol) / len(mol.atoms)


def n_halo(mol):
    # TODO: memoize
    """ The number of halogens.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """
    return sum(n_atom(mol, s) for s in ('F', 'Cl', 'Br', 'I'))


@requires_h_filled
def fract_halo(mol):
    """ The fraction of halogens.

     Args:
         mol (skchem.Mol):
             The molecule for which to calculate the descriptor.

     Returns:
         float
     """

    return n_halo(mol) / len(mol.atoms)


def n_hetero(mol):

    """ The number of heteroatoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return rdMolDescriptors.CalcNumHeteroatoms(mol)


def n_heavy(mol):

    """ The number of heavy atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return mol.GetNumHeavyAtoms()


def n_terminal(mol):

    """ The number of heavy atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return mol.atoms.is_terminal.sum()


@requires_h_filled
def n_bonds(mol):

    """ The number of bonds.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return len(mol.bonds)


@requires_h_depleted
def n_bonds_non_h(mol):

    """ The number of bonds between atoms other than hydrogen.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return len(mol.bonds)


@requires_h_depleted
def n_bonds_multiple(mol):

    """ The number of multiple bonds.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return (mol.bonds.order > 1).sum()


@requires_h_depleted
def sum_of_conventional_bond_orders(mol):

    """ The sum of conventional bond orders (h-depleted).

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int
    """

    return mol.bonds.order.sum()


def n_rotatable_bonds(mol):

    """ The number of rotatable bonds.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        int

    """
    return rdMolDescriptors.CalcNumRotatableBonds(mol)


@requires_h_depleted
def fract_rotatable_bonds(mol):
    # TODO: memoize
    """ The fraction of rotatable bonds.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return n_rotatable_bonds(mol) / n_bonds(mol)


@requires_h_depleted
def n_bond_order(mol, order=1):

    """ The number of bonds of order *i*.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        order (int):
            The order of bonds.
    Returns:
        int
    """
    return (mol.bonds.order == order).sum()


def fract_c_hybrid(mol, h_state='SP3'):

    """ The fraction of carbons that are in a certain hybridization state.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        h_state (str):
            The hybridization state for which to measure the fraction.
    Returns:
        float
    """

    carbs = mol.atoms.atomic_number == 6
    return (mol.atoms.hybridization_state == h_state)[carbs].sum() / carbs.sum()


def n_disconnected(mol):

    """ The number of disconnected fragments in the mol.

        Args:
            mol (skchem.Mol):
                The molecule for which to calculate the descriptor.

        Returns:
            int
        """

    return len(rdmolops.GetMolFrags(mol))


def total_charge(mol):

    """ The total charge of the molecule.

        Args:
            mol (skchem.Mol):
                The molecule for which to calculate the descriptor.

        Returns:
            float
        """

    return mol.atoms.formal_charge.sum()


def n_hba(mol):

    """ The number of h bond acceptors.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return rdMolDescriptors.CalcNumHBA(mol)


def n_hba(mol):

    """ The number of h bond donors.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return rdMolDescriptors.CalcNumHBD(mol)


def n_radical_electrons(mol):

    """ The number of radical electrons.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return Descriptors.NumRadicalElectrons(mol)


def n_valence_electrons(mol):

    """ The number of valence electrons.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return Descriptors.NumValenceElectrons(mol)


def heavy_atom_mol_wt(x):

    """ The molecular weight of only heavy atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return rdMolDescriptors.CalcExactMolWt(x, True)


def n_hbd_lipinski(x):

    """ The number of hydrogen bond donors according to Lipinski."""

    return rdMolDescriptors.CalcNumLipinskiHBD(x)


def n_hba_lipinski(x):

    """ The number of hydrogen bond acceptors according to Lipinski."""

    return rdMolDescriptors.CalcNumLipinskiHBA(x)


@requires_h_depleted
def n_paths(mol, length=1):

    """ The number of paths of length *l*. """

    return len(rdmolops.FindAllPathsOfLengthN(mol, length))


DESCRIPTORS = OrderedDict((
    ('mol_wt', molecular_weight),
    ('avg_mol_wt', average_molecular_weight),
    ('sum_vdw_vol', sum_van_der_waals_volume),
    ('sum_eneg', sum_electronegativity),
    ('sum_pol', sum_polarisability),
    ('sum_ion_energy', sum_ionisation_energy),
    ('mean_vdw_vol', mean_van_der_waals_volume),
    ('mean_eneg', mean_electronegativity),
    ('mean_pol', mean_polarisability),
    ('mean_ion_energy', mean_ionisation_energy),
    ('graph_density', graph_density),
    ('n_atoms', n_atoms),
    ('n_term', n_terminal),
    ('n_bonds', n_bonds),
    ('n_bonds_non_h', n_bonds_non_h),
    ('n_bonds_mult', n_bonds_multiple),
    ('sum_bond_order', sum_of_conventional_bond_orders),
    ('n_rot_bonds', n_rotatable_bonds),
    ('fract_rot_bonds', fract_rotatable_bonds),
    ('n_hyd', n_hyd),
    ('n_halo', n_halo),
    ('n_heavy', n_heavy),
    ('n_hetero', n_hetero),
    ('n_hba', n_hba),
    ('n_atoms', n_atoms),
    ('fract_halo', fract_halo),
    ('n_disconn', n_disconnected),
    ('total_charge', total_charge),
    ('n_rad', n_radical_electrons),
    ('n_val', n_valence_electrons),
    ('heavy_mol_wt', heavy_atom_mol_wt),
    ('n_hbd_lip', n_hbd_lipinski),
    ('n_hba_lib', n_hba_lipinski)
))

SYMBOLS = ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B')
HYBRIDS = ('SP3', 'SP2', 'SP')

DESCRIPTORS.update(('fract_c_{}'.format(hybrid), partial(fract_c_hybrid,
                                                         h_state=hybrid))
                   for hybrid in HYBRIDS)

DESCRIPTORS.update((('n_{}'.format(symbol), partial(n_atom, symbol=symbol))
                    for symbol in SYMBOLS))

DESCRIPTORS.update((('n_{}'.format(symbol), partial(fract_atom, symbol=symbol))
                    for symbol in ('H', 'C', 'N', 'O')))

DESCRIPTORS.update((('n_bond_{}'.format(order), partial(n_bond_order,
                                                       order=order))
                    for order in (1, 1.5, 2, 3)))

DESCRIPTORS.update((('n_paths_{}'.format(length), partial(n_paths, length=length))
                    for length in range(1, 7)))
