#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pytest
import json

from ...core import Mol

@pytest.fixture(name='m')
def ethane():
    return Mol.from_smiles('CC')


@pytest.fixture(name='m_c')
def ethane_copy():
    return Mol.from_smiles('CC')


@pytest.fixture(name='m2')
def propane():
    return Mol.from_smiles('CCC')


@pytest.fixture(name='m_b')
def bad_mol():
    return Mol.from_smarts('*~*')


@pytest.fixture(name='m_f')
def full_mol(m):
    return m.add_hs()


@pytest.fixture(name='d')
def ethane_dict():
    return {'m': [{'a': [{'l': 'C', 'x': -0.75, 'y': 0.0, 'z': 0.0},
                        {'l': 'C', 'x': 0.75, 'y': -0.0, 'z': 0.0}],
                  'b': [{'b': 0, 'e': 1, 'o': 1.0}]}]}


def test_smiles(m):
    assert len(m.atoms) == 2
    for a in m.atoms:
        assert a.symbol == 'C'


def test_substructure(m, m2):
    assert m in m2
    assert m2 not in m


def test_equality(m, m_c):
    assert m == m_c


def test_not_equals(m, m2):
    assert m != m2


def test_not_equals_other(m):
    assert m != 'notamol'


def test_binary(m):
    assert m == Mol.from_binary(m.to_binary())


def test_formula(m):
    assert m.to_formula() == 'C2H6'


def test_bad_formula(m_b):
    with pytest.raises(ValueError):
        m_b.to_formula()


def test_add_hs_inplace(m):
    with pytest.raises(NotImplementedError):
        m.add_hs(inplace=True)


def test_add_hs(m):
    assert len(m.add_hs().atoms) == 8


def test_remove_hs(m_f):
    assert len(m_f.remove_hs().atoms) == 2


def test_remove_hs_inplace(m):
    with pytest.raises(NotImplementedError):
        m.remove_hs(inplace=True)


def test_to_dict(m, d):
    assert m.to_dict() == d


def test_to_dict_not_conf(m, d):
    with pytest.raises(IndexError):
        m.to_dict(conformer_id=2)


def test_to_dict_other(m):
    with pytest.raises(NotImplementedError):
        m.to_dict(kind='doesntexist')


def test_to_json(m, d):
    assert json.loads(m.to_json()) == d


def test_to_json_other(m):
    with pytest.raises(NotImplementedError):
        m.to_json(kind='doesntexist')


def test_inchi_key(m):
    assert m.to_inchi_key() == 'OTMSDBZUPAUEDD-UHFFFAOYSA-N'


def test_inchi_key_fails(m_b):
    with pytest.raises(RuntimeError):
        m_b.to_inchi_key()


def test_inchi_key_not_avail(m):
    from rdkit.Chem import inchi

    inchi.INCHI_AVAILABLE, temp = False, inchi.INCHI_AVAILABLE
    with pytest.raises(ImportError):
        m.to_inchi_key()

    inchi.INCHI_AVAILABLE = temp


def test_copy(m):
    m_copy = m.copy()
    assert m_copy is not m
    assert m_copy == m


def test_contains_miss(m):
    with pytest.raises(NotImplementedError):
        'test' in m


def test_repr(m):
    assert repr(m) == '<Mol name="None" formula="C2H6" at {}>'.format(
        hex(id(m)))


def test_repr_smarts(m_b):
    assert repr(m_b) == '<Mol name="None" formula="unknown" at {}>'.format(
        hex(id(m_b)))



