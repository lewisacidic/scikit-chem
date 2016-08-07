#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.test

Tests for scikit-chem
"""

import pytest

chemaxon = pytest.mark.skipif(
    not pytest.config.getoption("--chemaxon"),
    reason="no chemaxon provided."
)