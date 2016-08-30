#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.test.test_features.caching

Tests for the caching module.
"""

import pytest

from ...features.descriptors import caching
from ..test_core import example_mol  # for m fixture


@pytest.fixture(name='cache')
def cache_fixture():
    return caching.Cache()


@pytest.fixture(name='func')
def func_fixture(cache):

    @cache
    def test_func(mol):
        return len(mol.atoms)

    return test_func


@pytest.fixture(name='kw_func')
def keyword_func_fixture(cache):

    @cache
    def test_func(mol, should_add_one=True):
        return len(mol.atoms) + should_add_one

    return test_func

def test_caches_func(m, cache, func):
    assert len(cache.cached) == 1

    @cache
    def second(mol):
        return 1

    assert len(cache.cached) == 2


def test_normal_use(m, cache, func):
    assert func(m) == 7


def test_caches_val(m, func):

    func(m)
    assert m.cache['test_func'][()] == 7


def test_uses_cached_val(m, cache, func):

    assert func(m) == 7

    # set a new value where it should be cached.  If we get this, the cached
    # value is being given
    m.cache['test_func'][()] = 100

    assert func(m) == 100


def test_dependency(m, cache):

    @cache
    def test_func(mol):
        return len(mol.atoms)

    @cache.inject(test_func)
    def test_dependency(mol, t_val):
        return t_val * 2

    assert test_dependency(m) == 14
    assert m.cache


def test_kwarg(m, cache, kw_func):

    @cache.inject(kw_func)
    def test_kw_dependency(mol, t_val, should_add_one=True):
        return t_val * 2

    assert test_kw_dependency(m) == 16
    assert test_kw_dependency(m, should_add_one=False) == 14

    assert cache