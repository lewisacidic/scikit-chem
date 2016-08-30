#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.features.descriptors.autocorrelation

Autocorrelation descriptors for scikit-chem.
"""

from functools import partial

from .caching import cache
from .fundamentals import adjacency_matrix, distance_matrix, atom_props

import numpy as np


@cache.inject(distance_matrix, atom_props)
def moreau_broto_autocorrelation(mol, dist_mat, prop, prop_name='atomic_mass',
                                 c_scaled=False, centred=False,
                                 ks=range(1, 9)):

    """ The Moreau-Broto autocorrelation.

    $$ ATS_k = \frac{1}{2} \hdot \sum_{i=1}^A \sum_{j=1}^A $$

    With special case $$ ATS_0 = \sum_{i=1}^A w_i^2.

    Where $A$ is the number of atoms, and $w$ is an atomic property.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        prop (str):
            The atomic property.

        c_scaled (bool):
            Whether the properties should be scaled against sp3 carbon.

        centred (bool):
            Whether the descriptor should be divided by the number of
            contributions (avoids dependence on molecular size).

        ks (iterable):
            The lags to calculate the descriptor over.

    Returns:
        float

    Examples:
        >>> import skchem
        >>> m = skchem.Mol.from_smiles('CC(O)CCO')
        >>> moreau_broto_autocorrelation(m, centred=False,
        ...                              prop_name='atomic_mass',
        ...                              ks=range(4))  # doctest: +ELLIPSIS
        array([ 1088.99...,   817.12...,   865.02...,   528.59...])
    """

    div = np.array([(dist_mat == k).sum() for k in ks]) if centred else 1

    return np.array([(1 if k == 0 else 0.5) * prop.dot(dist_mat == k).dot(prop)
                     for k in ks]) / div


@cache.inject(distance_matrix, atom_props)
def moran_coefficient(mol, dist_mat, prop, prop_name='atomic_mass',
                      c_scaled=False, ks=range(1, 9)):

    """ Moran coefficient for lags ks.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        prop (str):
            The atomic property.

        c_scaled (bool):
            Whether the properties should be scaled against sp3 carbon.

        centered (bool):
            Whether the descriptor should be divided by the number of
            contributions (avoids dependence on molecular size.

        ks (iterable):
            The lags to calculate the descriptor over.

    Returns:
        float
    """

    prop = prop - prop.mean()

    res = []

    for k in ks:
        geodesic = dist_mat == k
        num = prop.dot(geodesic).dot(prop) / geodesic.sum()
        denom = (prop ** 2).sum() / len(prop)
        res.append(num / denom)

    return np.array(res)


@cache.inject(distance_matrix, atom_props)
def geary_coefficient(mol, dist_mat, prop, prop_name='atomic_mass',
                      c_scaled=False, ks=range(1, 9)):

    """ The geary coefficient for *ks* lags.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        prop (str):
            The atomic property.

        c_scaled (bool):
            Whether the properties should be scaled against sp3 carbon.

        centered (bool):
            Whether the descriptor should be divided by the number of
            contributions (avoids dependence on molecular size.

        ks (iterable):
            The lags to calculate the descriptor over.

    Returns:
        float

    """

    res = []
    for k in ks:
        geodesic = dist_mat == k
        num = 0.5 * ((prop - prop[:, np.newaxis]) ** 2 * geodesic).sum() / geodesic.sum()
        denom = ((prop - prop.mean()) ** 2).sum() / (len(prop) - 1)
        res.append(num / denom)

    return np.array(res)


@cache
@cache.inject(adjacency_matrix, distance_matrix)
def galvez_matrix(mol, dist_mat, adj_mat):

    """ The galvez matrix.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the matrix.

    Returns:
        np.array
    """

    temp = dist_mat ** -2
    np.fill_diagonal(temp, 0)
    galvez_mat = temp.dot(adj_mat) + np.diag(mol.atoms.valence_vertex_degree)
    return galvez_mat


@cache
@cache.inject(galvez_matrix)
def charge_matrix(mol, galvez_mat):

    """ The charge matrix.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the matrix.

    Returns:
        np.array

    """

    ct_mat = galvez_mat - galvez_mat.T
    ct_mat[np.diag_indices_from(ct_mat)] = mol.atoms.depleted_degree

    return ct_mat


@cache
@cache.inject(charge_matrix, distance_matrix)
def topological_charge_index(mol, c_mat, dist_mat, ks=range(11)):

    """ The Galvez tologogical charge index for lags ks.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        ks (iterable):
            The lags for which to calculate the descriptor.

    Returns:
        np.array

    """
    return np.array([0.5 * np.abs(c_mat)[dist_mat == k].sum() for i in ks])


@cache.inject(topological_charge_index)
def mean_topological_charge_index(mol, tci, ks=range(11)):

    """ Mean topological charge index for lags ks.

    Args:
        mol (skchem.Mol):
            The molecule for which to calculate the descriptor.

        ks (iterable):
            The lags for which to calculate the descriptor.

    Returns:
        np.array
    """

    return tci / (len(mol.atoms) - 1)


@cache.inject(topological_charge_index)
def total_charge_index(mol, tci, ks=range(11)):

    return tci.sum() / (len(mol.atoms) - 1)


PROPS = ['atomic_mass', 'van_der_waals_volume', 'sanderson_electronegativity',
         'polarisability', 'ionisation_energy', 'intrinsic_state']

KS = range(1, 9)

DESCRIPTORS = [partial(moreau_broto_autocorrelation, k=k, prop=p, centered=c)
               for k in KS for c in (False, True) for p in PROPS]

FS = (moran_coefficient, geary_coefficient)

DESCRIPTORS += [partial(f, k=k, prop=p) for f in FS for k in KS for p in PROPS]

__all__ = ['moreau_broto_autocorrelation', 'moran_coefficient',
           'geary_coefficient', 'topological_charge_index',
           'mean_topological_charge_index', 'total_charge_index']
