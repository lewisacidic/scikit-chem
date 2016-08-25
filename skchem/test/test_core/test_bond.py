#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pytest
import numpy as np
import pandas as pd

from ...core import bond, Mol
from . import example_mol  #provides 'm' fixture


@pytest.fixture(name='b')
def example_bond(m):
    return m.bonds[0]


@pytest.fixture(name='bwp')
def example_bond_with_props(b):
    b.props['test'] = 'value'
    return b


def test_len(m):
    assert len(m.bonds) == 6


def test_out_of_range(m):
    with pytest.raises(IndexError):
        m.bonds[100]


def test_reverse_index(m):
    assert m.bonds[-1].order == 1

def test_slice(m):
    assert len(m.bonds[[1, 4]]) == 2


def test_repr(b):
    assert repr(b) ==  '<Bond type="O-C" at {}>'.format(hex(id(b)))


def test_owner(m):
    # rdkit gives a copy of the object, so cant test for identity
    assert m.bonds[0].owner.to_smiles() == m.to_smiles()


def test_to_dict(b):
    assert b.to_dict() == {'b': 0, 'e':1, 'o': 1}


def test_index(m):
    assert m.bonds.index.equals(pd.RangeIndex(6, name='bond_idx'))


def test_all_params_on_view():

    params = list(bond.Bond.__dict__.keys())

    for param in ('__doc__', '__repr__', '__str__', '__module__', 'atoms',
                  'props', 'owner', 'draw', 'to_dict'):

        params.remove(param)

    for param in params:
        assert hasattr(bond.BondView, param)


def test_atoms(b):
    assert len(b.atoms) == 2


def test_atom_idxs(b):
    assert b.atom_idxs == (0, 1)

test_data = [
    ('order', [1, 2, 1, 1, 1, 1]),
    ('stereo_symbol', ['NONE', 'NONE', 'NONE', 'NONE', 'NONE', 'NONE']),
    ('is_in_ring', [False, False, False, False, False, False]),
    ('is_conjugated', [True, True, False, False, False, False]),
    ('is_aromatic', [False, False, False, False, False, False])
]

params = pytest.mark.parametrize('param, expected', test_data)


@pytest.fixture(name='m_a')
def exampole_aromatic_mol():
    return Mol.from_smiles('c1ccNc1')


arom_test_data = [
    ('order', [1.5, 1.5, 1.5, 1.5, 1.5]),
    ('stereo_symbol', ['NONE', 'NONE', 'NONE', 'NONE', 'NONE']),
    ('is_in_ring', [True, True, True, True, True]),
    ('is_conjugated', [True, True, True, True, True]),
    ('is_aromatic', [True, True, True, True, True])
]

arom_params = pytest.mark.parametrize('param, expected', arom_test_data)


@params
def test_params_on_bond_view(m, param, expected):
    assert np.array_equal(getattr(m.bonds, param), expected)

@arom_params
def test_arom_params(m_a, param, expected):
    assert np.array_equal(getattr(m_a.bonds, param), expected)


@params
def test_params_on_bonds(m, param, expected):
    res = np.array([getattr(b, param) for b in m.bonds])
    assert np.array_equal(res, expected)


def test_props_keys_empty(b):
    assert len(b.props.keys()) == 0


def test_props_len_empty(b):
    assert len(b.props) == 0


def test_props_keys_full(bwp):
    assert len(bwp.props.keys()) == 1
    assert bwp.props.keys()[0] == 'test'


def test_props_len_full(bwp):
    assert len(bwp.props) == 1


def test_atom_idx_view(m):
    assert m.bonds.atom_idxs.shape == (len(m.bonds), 2)