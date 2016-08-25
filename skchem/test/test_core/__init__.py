#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.test.test_core

Tests for scikit-chem core functionality.
"""

import pytest
from rdkit.Chem.AllChem import Compute2DCoords, EmbedMultipleConfs

from ...core import Mol


@pytest.fixture(name='m')
def example_mol():
    m = Mol.from_smiles('[O-]C(=O)[C@H](F)CCl')
    Compute2DCoords(m, clearConfs=False)
    EmbedMultipleConfs(m, numConfs=10, clearConfs=False)
    return m
