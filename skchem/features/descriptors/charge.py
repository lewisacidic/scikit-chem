#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.features.descriptors.charge

Charge descriptors for scikit-chem.
"""

from collections import OrderedDict

import numpy as np
from rdkit.Chem import rdPartialCharges

from .caching import cache
from .fundamentals import geometric_matrix, adjacency_matrix

# TODO: consider using population analysis for 'proper' charge indices


def hstackable(res):
    return res if isinstance(res, np.ndarray) else np.array([res])


@cache
def gasteiger_charges(mol):

    """ The gasteiger partial charges for the atoms of the molecule.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptors.

    Returns:
        float
    """

    rdPartialCharges.ComputeGasteigerCharges(mol)
    return mol.atoms.props.pop('_GasteigerCharge')


@cache.inject(gasteiger_charges)
def max_partial_charge(mol, g_charges):

    """ The maximum gasteiger partial charge of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return g_charges.max()


@cache
@cache.inject(gasteiger_charges)
def min_partial_charge(mol, g_charges):

    """ The minimum gasteiger partial charge of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return g_charges.min()


@cache.inject(gasteiger_charges)
def total_positive_charge(mol, g_charges):

    """ The total positive gasteiger partial charges of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """


    return g_charges[g_charges > 0].sum()


@cache.inject(gasteiger_charges)
def total_negative_charge(mol, g_charges):

    """ The total negative gasteiger partial charges of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return g_charges[g_charges > 0].sum()


@cache.inject(gasteiger_charges)
def total_absolute_charge(mol, g_charges):

    """ The total absolute gasteiger partial charges of all atoms.

    This may be considered the electronic charge index (ECI).

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return np.abs(g_charges).sum()


@cache.inject(gasteiger_charges)
def mean_absolute_charge(mol, g_charges):

    """ The mean absolute gasteiger partial charges of all atoms.

    This may be considered the 'charge polarization'.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return np.abs(g_charges).mean()


@cache.inject(gasteiger_charges)
def total_squared_charge(mol, g_charges):

    """ The total of the squared  gasteiger partial charges of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return np.power(g_charges, 2).sum()


@cache.inject(gasteiger_charges, adjacency_matrix)
def submolecular_polarity_parameter(mol, g_charges, a_mat):

    """ The submolecular polarity parameter.

    This is the greatest difference in partial charges of bonded atoms in the
    molecule.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return (np.abs(g_charges - g_charges[:, np.newaxis]) * a_mat).max()


@cache
@cache.inject(gasteiger_charges, geometric_matrix)
def _topographic_electron_matrix(mol, g_charges, g_mat, conformer=-1):

    """ The matrix of topographic charge interactions.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    g_mat = g_mat ** -2
    np.fill_diagonal(g_mat, 0)
    return np.abs(g_charges - g_charges[:, np.newaxis]) * g_mat


@cache
@cache.inject(_topographic_electron_matrix)
def topographic_electronic_descriptor(mol, tem):

    """ The topolographic electronic descriptor.

    The sum of the differences in charges of all atoms in the molecule, scaled
    by the inverse square of their distances, i.e. a  measure of the coulombic
    interactions.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """
    return 0.5 * tem.sum()


@cache
@cache.inject(_topographic_electron_matrix, adjacency_matrix)
def br_topographic_electronic_descriptor(mol, tem, a_mat):

    """ The bond restricted topographic electronic descriptor.

    As the topographic electronic descriptor, except restricted only to bonded
    atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float

    See Also:
        skchem.features.descriptors.charge.topographic_electronic_descriptor
    """


    return 0.5 * (tem * a_mat).sum()


@cache.inject(topographic_electronic_descriptor, min_partial_charge)
def pcw_topological_electronic_index(mol, tem, q_min):

    """ Partial charge weighted topological electronic index.

    The topographic electronic index, scaled by the maximum negative partial
    charge in the molecule.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float

    See Also:
        skchem.features.descriptors.charge.topographic_electronic_descriptor
    """

    return np.abs(tem / q_min)


@cache.inject(br_topographic_electronic_descriptor, min_partial_charge)
def br_pcw_topological_electronic_index(mol, brtem, q_min):

    """ Bond restricted partial charge weighted topological electronic index.

    The partial charge weighted topological electronic index, with the
    summation restricted to bonded atom pairs.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float

    See Also:
        skchem.features.descriptors.charge.topographic_electronic_descriptor
    """
    return np.abs(brtem / q_min)


@cache.inject(gasteiger_charges, adjacency_matrix)
def local_dipole_index(mol, g_charges, a_mat):

    # cancel out 0.5s
    return (a_mat * np.abs(g_charges - g_charges[:, np.newaxis])).sum() / a_mat.sum()


@cache.inject(gasteiger_charges)
def min_abs_partial_charge(mol, g_charges):

    """ The minimum absolute gasteiger partial charge of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return np.abs(g_charges).min()


@cache.inject(gasteiger_charges)
def max_abs_partial_charge(mol, g_charges):

    """ The maximum absolute gasteiger partial charge of all atoms.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

    Returns:
        float
    """

    return np.abs(g_charges).max()


all_feats = OrderedDict(
    (('max_partial_charge', max_partial_charge),
     ('min_partial_charge', min_partial_charge),
     ('total_positive_charge', total_positive_charge),
     ('total_negative_charge', total_negative_charge),
     ('total_absolute_charge', total_absolute_charge),
     ('total_squared_charge', total_squared_charge),
     ('submolecular_polarity_parameter', submolecular_polarity_parameter),
     ('topographic_electronic_descriptor', topographic_electronic_descriptor),
     ('br_topographic_electronic_descriptor',
      br_topographic_electronic_descriptor),
     ('pcw_topological_electronic_index', pcw_topological_electronic_index),
     ('br_pcw_topological_electronic_index',
      br_pcw_topological_electronic_index),
     ('local_dipole_index', local_dipole_index),
     ('min_abs_partial_charge', min_abs_partial_charge),
     ('max_abs_partial_charge', max_abs_partial_charge))
)


__all__ = ['gasteiger_charges', 'max_partial_charge', 'min_partial_charge',
           'total_positive_charge', 'total_negative_charge',
           'total_absolute_charge', 'total_squared_charge',
           'submolecular_polarity_parameter',
           'topographic_electronic_descriptor',
           'br_topographic_electronic_descriptor',
           'pcw_topological_electronic_index',
           'br_pcw_topological_electronic_index',
           'local_dipole_index',
           'min_abs_partial_charge', 'max_abs_partial_charge']