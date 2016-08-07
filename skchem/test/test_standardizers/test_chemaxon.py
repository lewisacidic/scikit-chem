#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pandas as pd
import pytest

from ... import standardizers
from ...core import Mol
from .. import chemaxon

@pytest.fixture
def s():
    return standardizers.ChemAxonStandardizer()

@pytest.fixture
def m():
    return Mol.from_smiles('CCC.CC')

@chemaxon
def test_on_mol(s, m):
    m_s = s.transform(m)
    assert len(m.atoms) > len(m_s.atoms)

@chemaxon
def test_on_series(s, m):
    ser = pd.Series([m], index=['mol1'], name='structure')
    ser_s = s.transform(ser)
    assert len(ser_s) == len(ser)
    assert len(ser['mol1'].atoms) > len(ser_s['mol1'].atoms)
