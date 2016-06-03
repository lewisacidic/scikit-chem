#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import skchem
import pytest

@pytest.fixture
def m():
    return skchem.Mol.from_smiles('CC')
@pytest.fixture
def m_c():
    return skchem.Mol.from_smiles('CC')
@pytest.fixture
def m2():
    return skchem.Mol.from_smiles('CCC')

def test_smiles(m):
    assert len(m.atoms) == 2
    for a in m.atoms:
        assert a.element == 'C'

def test_substructure(m, m2):
    assert m in m2
    assert m2 not in m

def test_equality(m, m_c):
    assert m == m_c

def test_binary(m):
    assert m == skchem.Mol.from_binary(m.to_binary())
