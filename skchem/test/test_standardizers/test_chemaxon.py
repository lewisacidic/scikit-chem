#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import numpy as np
import pandas as pd
import pytest

import skchem

@pytest.fixture
def s():
    return skchem.standardizers.ChemAxonStandardizer()

@pytest.fixture
def m():
    return skchem.Mol.from_smiles('CCC.CC')

def test_on_mol(s, m):
    m_s = s.transform(m)
    assert len(m.atoms) > len(m_s.atoms)

def test_on_series(s, m):
    ser = pd.Series([m], index=['mol1'], name='structure')
    ser_s = s.transform(ser)
    assert len(ser_s) == len(ser)
    assert len(ser['mol1'].atoms) > len(ser_s['mol1'].atoms)
