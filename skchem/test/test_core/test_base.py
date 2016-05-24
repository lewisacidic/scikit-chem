#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pytest

import skchem

@pytest.fixture
def m():
    return skchem.Mol.from_smiles('CC')

@pytest.fixture
def mk(m):
    m.props['test'] = 'value'
    return m

@pytest.fixture
def ma(m):
    m.atoms.props['test'] = ['spam', 'eggs']
    return m

@pytest.fixture
def mkl(mk):
    mk.props['testerino']= 'valuerino'
    mk.props['testz'] = 'valz'
    return mk

@pytest.fixture
def a():
    return skchem.core.Atom('C')


def test_set_prop(mk):
    assert mk.props['test'] == 'value'

def test_miss(mk):
    with pytest.raises(KeyError):
        mk.props['notest']

def test_keys(mk):
    assert 'test' in mk.props.keys()

def test_len(mk):
    assert len(mk.props) == 1

def test_len2(mkl):
    assert len(mkl.props) == 3

def test_items(mk):
    assert ('test', 'value') in mk.props.items()

def test_get(mk):
    assert mk.props.get('test') == 'value'
    assert mk.props.get('notest') is None
    assert mk.props.get('notest', 'ham') == 'ham'

def test_pop(mk):
    assert mk.props.pop('test') == 'value'
    assert len(mk.props) == 0

def test_pop2(mkl):
    assert mkl.props.pop('test') == 'value'
    assert len(mkl.props) == 2

def test_remove(mk):
    mk.props.remove('test')
    assert len(mk.props) == 0

def test_remove2(mkl):
    mkl.props.remove('test')
    assert len(mkl.props) == 2

def test_clear(mk):
    mk.props.clear()
    assert len(mk.props) == 0

def test_clear2(mkl):
    mkl.props.clear()
    assert len(mkl.props) == 0

def test_delete(mkl):
    del mkl.props['test']
    assert len(mkl.props) == 2

def test_int_prop(m):
    m.props['test'] = 1
    assert m.props['test'] == 1

def test_float_prop(m):
    m.props['test'] = 3.142
    assert m.props['test'] == 3.142

def test_warns_key_not_string(m):
    with pytest.warns(UserWarning):
        m.props[4] = 'test'
    assert m.props['4'] == 'test'

def test_atom_prop(a):
    a.props['test'] = 'value'
    assert a.props['test'] == 'value'

def test_av_prop(ma):
    assert ma.atoms.props['test'] == ['spam', 'eggs']

def test_av_miss(ma):
    with pytest.raises(KeyError):
        ma.atoms.props['notest']

def test_av_keys(ma):
    assert 'test' in ma.atoms.props.keys()

def test_av_items(ma):
    assert ('test', ['spam', 'eggs']) in ma.atoms.props.items()

def test_av_len(ma):
    assert len(ma.atoms.props) == 1

def test_av_del(ma):
    del ma.atoms.props['test']
    assert len(ma.atoms.props) == 0
