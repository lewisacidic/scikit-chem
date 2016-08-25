#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.test

Tests for scikit-chem
"""

import pytest

class FakeConfig(object):
    def getoption(self, arg):
        pass

if not hasattr(pytest, 'config'):
    pytest.config = FakeConfig()

with_chemaxon = pytest.mark.skipif(
    not pytest.config.getoption('--with-chemaxon'),
    reason='no chemaxon provided.')

slow = pytest.mark.skipif(
    not pytest.config.getoption('--run-slow'),
    reason='test runs slowly.')