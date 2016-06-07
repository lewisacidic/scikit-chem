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

def test_on_df(s, m):
    df = pd.DataFrame([[m, 'spam', 'eggs']],
                      index=['mol1'],
                      columns=['structure', 'prop1', 'prop2'])
    df_s = s.transform(df)
    assert np.array_equal(df.index, df_s.index)
    assert np.array_equal(df.columns, df_s.columns)
    assert len(df.ix['mol1', 'structure'].atoms) > len(df_s.ix['mol1', 'structure'].atoms)
