#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pytest

from rdkit.Chem.rdDepictor import Compute2DCoords
import numpy as np

from . import example_mol  #for m fixture
from ...core import Mol
from ...core.conformer import Conformer


@pytest.fixture(name='c')
def two_d_conformer(m):
    return m.conformers[0]


@pytest.fixture(name='c3d')
def three_d_conformer(m):
    return m.conformers[1]


@pytest.fixture(name='l')
def linear_mol():
    l = Mol.from_smiles('O=C=O')
    Compute2DCoords(l)
    return l


@pytest.fixture(name='rot_mat')
def random_rotation_mat():
    rand = np.random.uniform(1, np.pi)
    return np.array([[1, 0, 0],
                     [0, np.cos(rand), -np.sin(rand)],
                     [0, np.sin(rand), np.cos(rand)]])


@pytest.fixture(name='c_l')
def linear_conformer(l):
    return l.conformers[0]


@pytest.fixture(name='c_rot')
def linear_randomly_rotated(c_l, rot_mat):
    rand = np.random.uniform(1, np.pi)
    c_l.positions.dot(rot_mat)
    return c_l


@pytest.fixture(name='c_pos')
def linear_randomly_translated(c_l):
    c_l.positions += 5 * np.random.randn(3)
    return c_l


@pytest.fixture(name='c_rand')
def linear_randomly_rotated_and_translated(c_l, rot_mat):
    c_l.positions += 5 * np.random.randn(3)
    c_l.positions = c_l.positions.dot(rot_mat)
    return c_l


@pytest.fixture(name='g')
def example_geom_matrix(m):
    return np.random.randn(3 * len(m.atoms)).reshape(len(m.atoms), 3)


# rdkit owner returns a copy rather that same object.
def test_owner(m, c):
    assert c.owner.to_smiles() == m.to_smiles()


def test_positions(m, c):
    assert np.array_equal(c.positions.shape, (len(m.atoms), 3))


def test_set_positions(c, g):
    c.positions = g
    assert np.array_equal(c.positions, g)


def test_set_2d_positions(c, g):
    c.positions = g[:, :2]

    assert np.array_equal(c.positions[:, :2], g[:, :2])
    assert (c.positions[:, 2] == 0).all()


def test_geometric_centre(c_pos):
    expected = c_pos.positions.mean(axis=0)
    assert np.array_equal(c_pos.geometric_centre, expected)


def test_centre_of_mass(l, c_pos):
    masses = l.atoms.atomic_mass
    expected = masses.dot(c_pos.positions) / masses.sum()
    assert np.array_equal(c_pos.centre_of_mass, expected)


def test_centre_representation(c_pos):
    c_pos.centre_representation(centre_of_mass=True)
    assert np.allclose(c_pos.centre_of_mass, (0, 0, 0))

    c_pos.centre_representation(centre_of_mass=False)
    assert np.allclose(c_pos.geometric_centre, (0, 0, 0))


def test_align(c_rot):
    c_rot.align_with_principal_axes()
    x, y, z = c_rot.positions.T
    assert np.allclose(y, np.zeros(len(y)))
    assert np.allclose(z, np.zeros(len(z)))


def test_canonicalize(c_rand):
    c_rand.canonicalize()
    x, y, z = c_rand.positions.T
    assert np.allclose(y, np.zeros(len(y)))
    assert np.allclose(z, np.zeros(len(z)))
    assert np.allclose(c_rand.centre_of_mass, (0, 0, 0))


def test_id(c):
    assert c.id == c.GetId()


def test_set_id(c):
    with pytest.warns(UserWarning):
        c.id = 1
    assert c.GetId() == 1


def test_is_3d(c, c3d):
    assert c.is_3d == False
    assert c3d.is_3d == True


def test_set_3d(c3d):
    assert c3d.is_3d == True
    c3d.is_3d = False
    assert c3d.is_3d == False
    assert all(c3d.positions[:, 2] == 0)


def test_repr(c):
    assert repr(c) == '<Conformer id="0" at {}>'.format(hex(id(c)))


# view tests
def test_get_conformer(m):
    assert isinstance(m.conformers[0], Conformer)


def test_miss(m):
    with pytest.raises(IndexError):
        m.conformers[11]


def test_add(m):
    c = Conformer(len(m.atoms))
    m.conformers[11] = c
    assert len(m.conformers) == 12
    assert c.id == 11


def test_remove(m):
    del m.conformers[0]
    assert len(m.conformers) == 10


def test_append(m):
    c = Conformer(len(m.atoms))
    m.conformers.append(c)
    assert len(m.conformers) == 12
    assert c.id == 11


def test_append_2d(m):
    m.conformers.append_2d()
    assert len(m.conformers) == 12
    assert m.conformers[11].is_3d == False


def test_append_3d(m):
    m.conformers.append_3d()
    assert len(m.conformers) == 12
    assert m.conformers[11].is_3d == True


def test_slicing(m):
    assert len(m.conformers[[0, 1]]) == 2
    assert len(m.conformers[:2]) == 2
    # cant check for equality as rdkit produces copies
    assert np.array_equal(m.conformers[-1].positions,
                          m.conformers[0].positions)


def test_positions(m):
    assert m.conformers.positions.shape == (len(m.conformers), len(m.atoms), 3)


def test_is_3d(m):
    assert np.array_equal(m.conformers.is_3d, [False] + [True] * 10)
