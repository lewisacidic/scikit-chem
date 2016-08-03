#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pytest
import pandas as pd
import numpy as np

from ...core import Mol
from ...filters import Filter, AtomNumberFilter, MassFilter, SMARTSFilter, ElementFilter

@pytest.fixture
def m():
    return Mol.from_smiles('CC', name='test')

@pytest.fixture
def ms():
    return [Mol.from_smiles(s, name=n) for s, n in \
            zip(('C', 'CC', 'CCC'), ('a', 'b', 'c'))]

@pytest.fixture
def f():
    return Filter(lambda mol: len(mol.atoms) > 1)

def test_takes_mol(m, f):
    assert f.transform(m) == True

def test_takes_mol_transform(m, f):
    assert f.transform(m) == True

def test_takes_list(m, f):
    assert np.array_equal(f.transform([m]), [True])

def test_takes_dict(m, f):
    assert np.array_equal(f.transform({'name': m}), [True])

def test_takes_ser(m, f):
    assert np.array_equal(f.transform(pd.Series({'name': m})), [True])

def test_filter(ms, f):
    assert len(f.filter(ms)) == 2
