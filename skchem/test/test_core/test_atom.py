#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


import pytest
import numpy as np
import pandas as pd

from ...core import atom, Mol


from . import example_mol


@pytest.fixture(name='plat')
def plat_mol():
    return Mol.from_smiles('[Pt]')


@pytest.fixture(name='a')
def example_atom(m):
    return m.atoms[0]


@pytest.fixture(name='dummy_props')
def dummy_properties():
    return ['a', 'b', 'c', 'd', 'e', 'f', 'g']

@pytest.fixture(name='mwp')
def mol_with_props(m, dummy_props):
    m.atoms.props['test'] = dummy_props
    return m


@pytest.fixture(name='awp')
def atom_with_props(a):
    a.props['test'] = 'value'
    return a

xfail = pytest.mark.xfail

# rdkit doesn't model COO atoms as delocalized.
test_data = [
    ('symbol', ['O', 'C', 'O', 'C', 'F', 'C', 'Cl']),
    ('hexcode', ['#FF0D0D', '#909090', '#FF0D0D', '#909090', '#90E050',
                 '#909090', '#1FF01F']),
    ('atomic_mass', [15.999, 12.011, 15.999, 12.011, 18.998, 12.011, 35.453]),
    ('n_total_hs', [0, 0, 0, 1, 0, 2, 0]),
    ('formal_charge', [-1, 0, 0, 0, 0, 0, 0]),
    ('n_hs', [0, 0, 0, 1, 0, 2, 0]),
    ('hybridization_state', ['SP2', 'SP2', 'SP2', 'SP3', 'SP3', 'SP3', 'SP3']),
    ('electron_affinity', [1.461, 1.596, 1.461, 1.596, 3.399, 1.596, 3.617]),
    ('n_lone_pairs', [3, 0, 2, 0, 3, 0, 3]),
    ('degree', [1, 3, 1, 3, 1, 2, 1]),
    ('mcgowan_parameter', [12.43, 16.35, 12.43, 16.35, 10.48, 16.35, 20.94]),
    ('kier_hall_alpha_contrib', [-0.195, -0.13, -0.195, 0, -0.065, 0, 0.286]),
    ('kier_hall_electronegativity', [1.25, 0.25, 1.25, 0., 1.5, 0., 0.667]),
    ('is_terminal', [True, False, True, False, True, False, True]),
    ('n_val_electrons', [6, 4, 6, 4, 7, 4, 7]),
    ('polarisability', [0.80, 1.76, 0.80, 1.76, 0.56, 1.76, 2.18]),
    ('principal_quantum_number', [2, 2, 2, 2, 2, 2, 3]),
    ('chiral_tag', [0, 0, 0, 2, 0, 0, 0]),
    ('n_pi_electrons', [0, 1, 1, 0, 0, 0, 0]),
    ('explicit_valence', [1, 4, 2, 4, 1, 2, 1]),
    ('full_degree', [1, 3, 1, 4, 1, 4, 1]),
    ('valence', [1, 4, 2, 4, 1, 4, 1]),
    ('n_instanced_hs', [0, 0, 0, 0, 0, 0, 0]),
    ('atomic_number', [8, 6, 8, 6, 9, 6, 17]),
    ('n_implicit_hs', [0, 0, 0, 1, 0, 0, 0]),
    ('depleted_degree', [1, 3, 1, 3, 1, 2, 1]),
    ('intrinsic_state', [7.,  1.667,  7.,  1.333,  8., 1.5,  4.111]),
    ('covalent_radius', [0.62, 0.67, 0.62, 0.77, 0.72, 0.77, 0.99]),
    ('van_der_waals_radius', [1.4, 1.75, 1.4, 1.75, 1.30, 1.75, 1.75]),
    ('n_explicit_hs', [0, 0, 0, 0, 0, 2, 0]),
    ('van_der_waals_volume', [11.494, 22.449, 11.494, 22.449, 9.2030, 22.449,
                              22.449]),
    ('is_aromatic', [False, False, False, False, False, False, False]),
    ('sanderson_electronegativity', [3.65, 2.75, 3.65, 2.75, 4.0, 2.75, 3.48]),
    ('implicit_valence', [0, 0, 0, 0, 0, 2, 0]),
    ('pauling_electronegativity', [3.44, 2.55, 3.44, 2.55, 3.98, 2.55, 3.16]),
    ('ionisation_energy', [13.61, 11.25, 13.61, 11.25, 17.41, 11.25, 13.01]),
    ('cahn_ingold_prelog', [None, None, None, 'S', None, None, None]),
    ('is_in_ring', [False, False, False, False, False, False, False])
]

atom_params = pytest.mark.parametrize('field, expected', test_data)

valence_test_data = [
    ('degree', [1, 3, 1, 2, 2, 1]),
    ('valence_degree', [1, 4, 5, 3, 2, 3]),
    xfail(('bond_vertex_degree', [1, 4, 1, 3, 2, 1])),
    xfail(('valence_state_indicator', [2, 7, 6, 5, 4, 4])),
    xfail(('perturbation_delta_value', [1.4, 4, 5, 4.9, 3.6, 2.6, 3.2])),
    xfail(('kupchik_vertex_degree', [1, 3, 5.274, 2, 2, 3.08])),
    xfail(('hu_xu_vertex_degree', [2.449, 7.348, 2.828, 4.898, 4.898, 2.646])),
    xfail(('alikanidi_vertex_degree',
           [2.723, 25.542, 2.723, 10.892, 11.247, 2.723])),
    ('intrinsic_state', [2, 1.667, 6, 2, 1.5, 4]),
    xfail(('ren_vertex_degree', [1, 3.2, 1.167, 2.25, 2, 1.25])),
    xfail(('li_vertex_degree', [1, 4, 7.5, 3, 2, 3.75])),
    xfail(('yang_vertex_degree', [0.588, 1.765, 1.302, 1.324, 1.176, 1.216])),
    xfail(('madan_chemical_degree',
           [1, 3.332, 1, 2, 2.167, 1])),
    xfail(('extended_madan_degree', [1, 4.332, 1, 3, 2.167, 1])),
    xfail(('ct_vertex_degree', [1, 3.46, 1.08, 2.378, 2.041, 1.041])),
    xfail(('z_delta_number', [2, 2, 3, 2, 2, 2.5]))
]

valence_params = pytest.mark.parametrize('field, expected', valence_test_data)


def test_all_params_on_view():

    params = list(atom.Atom.__dict__.keys())

    for param in ('__doc__', '__repr__', '__str__', '__module__', '_cov_dict',
                  'props', 'owner', 'bonds', 'neighbours'):
        params.remove(param)

    for param in params:
        assert hasattr(atom.AtomView, param)


# Atom tests

# some bad tests - rdkit creates copies of objects rather than the objects
# themselves, so this is the best we can do.
def test_neighbours(m):

    assert m.atoms[0].neighbours()[0].symbol == m.atoms[1].symbol


def test_bonds(m):

    bs = m.atoms[0].bonds
    assert len(bs) == 1
    assert bs[0].order == 1


def test_owner(m):

    assert m.atoms[0].owner.to_smiles() == m.to_smiles()


@atom_params
def test_params_by_atoms(m, field, expected):

    res = np.array([getattr(a, field) for a in m.atoms])

    if (res.dtype.type is np.str_) or (res.dtype.type is np.object_):
        assert np.array_equal(res, np.array(expected))
    else:
        assert np.allclose(res, expected, atol=0.01)


@valence_params
def test_valence_vertex_degree(val_m, field, expected):
    assert np.allclose([getattr(a, field) for a in val_m.atoms],
                       expected, atol=0.001)


def test_rare_covalent(plat):
    assert np.allclose(plat.atoms.covalent_radius, [1.28], atol=0.01)


def test_h_kier_hall(m):
    assert m.add_hs().atoms[-1].kier_hall_electronegativity == -0.2


def test_repr(a):
    assert repr(a) == '<Atom element="O" at {}>'.format(hex(id(a)))


def test_str(a):
    assert str(a) == 'O'


@pytest.fixture
def val_m():
    return Mol.from_smiles('CC(O)=CCN')


# test Atom props

def test_no_props_keys(a):
    assert len(a.props.keys()) == 0


def test_no_props_len(a):
    assert len(a.props) == 0


def test_raises_key_error(a):
    with pytest.raises(KeyError):
        a.props['test']


def test_keys(awp):
    assert len(awp.props.keys()) == 1


def test_getattr(awp):
    assert awp.props['test'] == 'value'


def test_get(awp):
    assert awp.props.get('test') == 'value'


def test_get_with_default(awp):
    assert awp.props.get('test', None) == 'value'


def test_get_miss_with_default(awp):
    assert awp.props.get('not_prop', 'other_val') == 'other_val'


def test_props_pop(awp):
    assert len(awp.props) == 1
    assert awp.props.pop('test') == 'value'
    assert len(awp.props) == 0


def test_props_pop_default(awp):
    assert awp.props.pop('test', 'default') == 'value'
    assert awp.props.pop('missing', 'default') == 'default'


def test_props_string(awp):
    assert str(awp.props) == '{\'test\': \'value\'}'

# AtomView tests

@atom_params
def test_params_by_atom_view(m, field, expected):

    res = getattr(m.atoms, field)

    if (res.dtype.type is np.str_) or (res.dtype.type is np.object_):
        assert np.array_equal(res, np.array(expected).astype(str))
    else:
        assert np.allclose(res, expected, atol=0.01)


@valence_params
def test_valence_vertex_degree_by_atom_view(val_m, field, expected):
    assert np.allclose(getattr(val_m.atoms, field), expected, atol=0.001)


def test_len(m):
    assert len(m.atoms) == 7


def test_slice(m):
    assert len(m.atoms[:4]) == 4
    assert [a.symbol for a in m.atoms[::2]] == ['O', 'O', 'F', 'Cl']
    assert [a.symbol for a in m.atoms[[1, 3, 5]]] == ['C', 'C', 'C']
    assert [a.symbol for a in m.atoms[(0, 2)]] == ['O', 'O']
    assert [a.symbol for a in m.atoms[
        (False, False, False, False, True, False, True)]] == ['F', 'Cl']


def test_index_too_high(m):
    with pytest.raises(IndexError):
        m.atoms[1000]

# props
def test_set_prop(m, dummy_props):
    m.atoms.props['test'] = dummy_props
    assert np.array_equal(m.atoms.props['test'], dummy_props)
    assert m.atoms[0].props['test'] == dummy_props[0]


def test_set_prop_dict(m):
    m.atoms.props['test'] = {1: 'value'}
    expected = [None, 'value'] + [None] * 5
    assert np.array_equal(m.atoms.props['test'], expected)


def test_view_props_keys_empty(m):
    assert len(m.atoms.props.keys()) == 0


def test_view_props_len_empty(m):
    assert len(m.atoms.props) == 0


def test_view_props_keys_full(mwp):
    assert len(mwp.atoms.props.keys()) == 1


def test_view_props_len_full(mwp):
    assert len(mwp.atoms.props) == 1


def test_view_props_pop(mwp, dummy_props):
    assert np.array_equal(mwp.atoms.props.pop('test'), dummy_props)
    with pytest.raises(KeyError):
        mwp.atoms.props.pop('test')


def test_view_props_iter(mwp):
    len([prop for prop in mwp.atoms.props]) == 1


def test_view_warns_reserved(m, dummy_props):
    with pytest.warns(UserWarning):
        m.atoms.props['_test'] = dummy_props


def test_view_warns_obj(m):
    with pytest.warns(UserWarning):
        m.atoms.props['test'] = [{}] * len(m.atoms)


def test_view_index(m):
    assert m.atoms.index.equals(pd.RangeIndex(7, name='atom_idx'))
