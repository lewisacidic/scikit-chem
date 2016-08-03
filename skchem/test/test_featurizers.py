
from ..descriptors import AtomFeaturizer
from ..core import Mol
import pytest

import numpy as np
import pandas as pd


@pytest.fixture
def af():
    return AtomFeaturizer()


@pytest.fixture
def m():
    return Mol.from_smiles('CCO')


@pytest.fixture
def a(m):
    return m.atoms[0]


@pytest.fixture
def s(m):
    return pd.Series([m, m], index=(1, 2))


def test_af(af):
    assert len(af.minor_axis) == 37


def test_on_a(af, a):
    feats = af.transform(a)
    assert np.array_equal(feats.shape, (len(af.minor_axis),))
    assert feats.index.name == 'atom_features'
    assert feats['is_C'] == True
    assert feats['is_SP2_hybridized'] == False


def test_on_m(af, m):
    feats = af.transform(m)
    assert np.array_equal(feats.shape, (len(m.atoms), len(af.minor_axis)))
    assert feats.index.name == 'atom_idx'
    assert feats.columns.name == 'atom_features'
    assert feats.ix[0, 'is_C'] == True
    assert feats.ix[2, 'is_C'] == False
    assert feats.ix[2, 'is_O'] == True


def test_on_ser(af, s):
    print(s)
    feats = af.transform(s)
    print(feats)
    assert np.array_equal(feats.shape, (len(s), af.max_atoms, len(af.minor_axis)))
    assert feats.ix[1, 0, 'is_C'] == True
    assert feats.ix[2, 2, 'is_O'] == True