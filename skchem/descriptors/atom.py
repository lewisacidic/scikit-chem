#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.atom

Module specifying atom based descriptor generators.
"""

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors, rdPartialCharges
from rdkit.Chem.rdchem import HybridizationType

import functools

from .. import core
from ..data import PERIODIC_TABLE
from ..filters import ORGANIC

def element(a):

    """ Return the element """

    return a.GetSymbol()

def is_element(a, symbol='C'):

    """ Is the atom of a given element """
    return element(a) == symbol

element_features = {'is_{}'.format(e): functools.partial(is_element, symbol=e) for e in ORGANIC}

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

    return a.mass

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

def gasteiger_charge(a, force_calc=False):

    """ Hacky way of getting gasteiger charge """

    res = a.props.get('_GasteigerCharge', None)
    if res and not force_calc:
        return float(res)
    else:
        idx = a.GetIdx()
        m = a.GetOwningMol()
        rdPartialCharges.ComputeGasteigerCharges(m)
        return float(a.props['_GasteigerCharge'])

def electronegativity(a):

    return PERIODIC_TABLE.loc[a.atomic_number, 'pauling_electronegativity']

def first_ionization(a):

    return PERIODIC_TABLE.loc[a.atomic_number, 'first_ionisation_energy']

def group(a):

    return PERIODIC_TABLE.loc[a.atomic_number, 'group']

def period(a):

    return PERIODIC_TABLE.loc[a.atomic_number, 'period']

def is_hybridized(a, hybrid_type=HybridizationType.SP3):

    """ Hybridized as type hybrid_type, default SP3 """

    return str(a.GetHybridization()) is hybrid_type

hybridization_features = {'is_' + n + '_hybridized': functools.partial(is_hybridized, hybrid_type=n) for n in HybridizationType.names}

atom_features = {
    'atomic_number': atomic_number,
    'atomic_mass': atomic_mass,
    'formal_charge': formal_charge,
    'gasteiger_charge': gasteiger_charge,
    'electronegativity': electronegativity,
    'first_ionisation': first_ionization,
    'group': group,
    'period': period,
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
}
atom_features.update(element_features)
atom_features.update(hybridization_features)

class AtomFeatureCalculator(object):

    def __init__(self, features='all', max_atoms=50):
        if features == 'all':
            features = atom_features
        self.features = pd.Series(features)
        self.max_atoms = max_atoms

    @property
    def feature_names(self):
        return self.features.index

    def transform(self, obj):
        if isinstance(obj, core.Atom):
            return self._transform_atom(obj)
        elif isinstance(obj, core.Mol):
            return self._transform_mol(obj)
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_atom(self, atom):
        return self.features.apply(lambda f: f(atom))

    def _transform_mol(self, mol):
        return pd.DataFrame(self.transform(a) for a in mol.atoms)

    def _transform_series(self, data):
        X = np.repeat(np.nan,
                      len(self.feature_names) * self.max_atoms * len(data))\
                      .reshape((len(data), self.max_atoms, len(self.feature_names)))
        for i, mol in enumerate(data):
            X[i, :len(mol.atoms), :] = self._transform_mol(mol)
        res = pd.Panel(X, items=data.index)
        res.minor_axis = self.features.index
        res.major_axis.name = 'atom_idx'
        return res


    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)

class GraphDistanceCalculator(object):

    def __init__(self, max_atoms=50):
        self.max_atoms = max_atoms

    def transform(self, obj):
        if isinstance(obj, core.Atom):
            return self._transform_atom(obj)
        elif isinstance(obj, core.Mol):
            return self._transform_mol(obj)
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_mol(self, mol):
        res = np.repeat(np.nan, self.max_atoms ** 2).reshape(self.max_atoms, self.max_atoms)
        res[:len(mol.atoms), :len(mol.atoms)] = Chem.GetDistanceMatrix(mol)
        return res

    def _transform_series(self, ser):
        res = pd.Panel(np.array([self.transform(mol) for mol in ser]), items=ser.index)
        res.major_axis.name = res.minor_axis.name = 'atom_idx'

    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)
