#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import HybridizationType
import skchem

import functools
from skchem.data import periodic_table

from .filters import organic

def is_element(a, element=0):

    """ Is the atom of a given element """
    return a.GetSymbol() == element

element_features = {'is_{}'.format(e): functools.partial(is_element, element=e) for e in organic}

def is_h_acceptor(a):

    """ Is an H acceptor? """

    m = a.GetOwningMol()
    idx = a.GetIdx()
    return idx in [i[0] for i in Lipinski._HAcceptors(m)]

def is_h_donor(a):

    """ Is an H donor? """

    m = a.GetOwningMol()
    idx = a.GetIdx()
    return idx in [i[0] for i in Lipinski._HDonors(m)]

def is_hetero(a):

    """ Is a heteroatom? """

    m = a.GetOwningMol()
    idx = a.GetIdx()
    return idx in [i[0] for i in Lipinski._Heteroatoms(m)]

def atomic_number(a):

    """ Atomic number of atom """

    return a.GetAtomicNum()

def atomic_mass(a):

    """ Atomic mass of atom """

    return a.GetMass()

def explicit_valence(a):

    """ Explicit valence of atom """
    return a.GetExplicitValence()

def implicit_valence(a):

    """ Implicit valence of atom """

    return a.GetImplicitValence()

def valence(a):

    """ returns the valence of the atom """

    return explicit_valence(a) + implicit_valence(a)

def formal_charge(a):

    """ Formal charge of atom """

    return a.GetFormalCharge()

def is_aromatic(a):

    """ Boolean if atom is aromatic"""

    return a.GetIsAromatic()

def num_implicit_hydrogens(a):

    """ Number of implicit hydrogens """

    return a.GetNumImplicitHs()

def num_explicit_hydrogens(a):

    """ Number of explicit hydrodgens """

    return a.GetNumExplicitHs()

def num_hydrogens(a):

    """ Number of hydrogens """

    return num_implicit_hydrogens(a) + num_explicit_hydrogens(a)

def is_in_ring(a):

    """ Whether the atom is in a ring """

    return a.IsInRing()

def crippen_log_p_contrib(a):

    """ Hacky way of getting logP contribution. """

    idx = a.GetIdx()
    m = a.GetOwningMol()
    return Crippen._GetAtomContribs(m)[idx][0]

def crippen_molar_refractivity_contrib(a):

    """ Hacky way of getting molar refractivity contribution. """

    idx = a.GetIdx()
    m = a.GetOwningMol()
    return Crippen._GetAtomContribs(m)[idx][1]

def tpsa_contrib(a):

    """ Hacky way of getting total polar surface area contribution. """

    idx = a.GetIdx()
    m = a.GetOwningMol()
    return rdMolDescriptors._CalcTPSAContribs(m)[idx]

def labute_asa_contrib(a):

    """ Hacky way of getting accessible surface area contribution. """

    idx = a.GetIdx()
    m = a.GetOwningMol()
    return rdMolDescriptors._CalcLabuteASAContribs(m)[0][idx]

def is_hybridized(a, hybrid_type=HybridizationType.SP3):

    """ Hybridized as type hybrid_type, default SP3 """

    return a.GetHybridization() is hybrid_type

pt = periodic_table()

def electronegativity(a):

    return a.GetSymbol()


hybridization_features = {'is_' + n + '_hybridized': functools.partial(is_hybridized, hybrid_type=n) for n in HybridizationType.names}

atom_features = {
    'atomic_number': atomic_number,
    'atomic_mass': atomic_mass,
    'formal_charge': formal_charge,
    'valence': valence,
    'is_aromatic': is_aromatic,
    'num_hydrogens': num_hydrogens,
    'is_in_ring': is_in_ring,
    'log_p_contrib': crippen_log_p_contrib,
    'molar_refractivity_contrib': crippen_molar_refractivity_contrib,
    'is_h_acceptor': is_h_acceptor,
    'is_h_donor': is_h_donor,
    'is_heteroatom': is_hetero,
    'total_polar_surface_area_contrib': tpsa_contrib,
    'total_labute_accessible_surface_area': labute_asa_contrib,
    **element_features,
    **hybridization_features
}

class AtomFeatureCalculator(object):

    def __init__(self, features=atom_features):
        self.features = pd.Series(features)

    @property
    def feature_names(self):
        return self.features.index

    def calculate(self, obj):
        if isinstance(obj, skchem.core.Atom):
            return self._calculate_atom(obj)
        elif isinstance(obj, skchem.core.Mol):
            return self._calculate_mol(obj)

    def _calculate_atom(self, atom):
        return self.features.apply(lambda f: f(atom))

    def _calculate_mol(self, mol):
        return pd.DataFrame(self(a) for a in mol.atoms)

    def __call__(self, *args, **kwargs):
        return self.calculate(*args, **kwargs)
