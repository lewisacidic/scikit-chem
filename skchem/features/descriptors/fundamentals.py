#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.features.descriptors.fundamentals

Fundamental representations for descriptors.
"""


from rdkit.Chem import rdmolops
import numpy as np

from .caching import cache
from ...core import Mol


@cache
def atom_props(mol, prop_name='unweighted', c_scaled=False):

    """ Atom based properties."""

    if prop_name == 'unweighted':
        return np.ones(len(mol.atoms))
    else:
        props = getattr(mol.atoms, prop_name)
        if c_scaled:
            props /= getattr(Mol.from_smiles('CC').atoms[0], prop_name)
        return props




@cache
def distance_matrix(mol):

    """ The topological distance matrix. """

    return rdmolops.GetDistanceMatrix(mol)


@cache
def adjacency_matrix(mol):

    """ The topological adjacency matrix. """

    return rdmolops.GetAdjacencyMatrix(mol)


@cache
def bond_order_adjacency_matrix(mol):

    """ The bond order scaled topological adjacency matrix. """
    return rdmolops.GetAdjacencyMatrix(mol, useBO=1)


# probably this anymore (props, degrees)
@cache
def degrees(mol):

    return mol.atoms.degree


@cache
def geometric_matrix(mol, conformer=-1):
    return rdmolops.Get3DDistanceMatrix(mol, confId=conformer)


@cache
def molecular_matrix(mol, conformer=-1):
    return mol.conformers[conformer].positions


__all__ = ['distance_matrix', 'adjacency_matrix',
           'bond_order_adjacency_matrix', 'degrees', 'geometric_matrix',
           'molecular_matrix']