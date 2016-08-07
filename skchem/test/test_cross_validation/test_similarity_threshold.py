#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.tests.test_cross_validation.test_similarity_threshold

Tests for similarity threshold dataset partitioning functionality.
"""

import pytest

from scipy.spatial.distance import cdist
import numpy as np

from ...data import Diversity
from ...cross_validation import SimThresholdSplit


@pytest.fixture
def x():
    return Diversity.read_frame('feats/X_morg')

@pytest.fixture
def cv(x):
    return SimThresholdSplit(x, fper=None, block_width=500, n_jobs=1)

def test_split(cv, x):
    train, test = cv.split((8, 2))
    assert (1 - cdist(x[train], x[test]) > cv.threshold_).sum() == 0
    assert np.allclose([train.sum()], [len(x) * 0.8], rtol=0.05)

def test_k_fold(cv, x):
    kfold = [fold for fold in cv.k_fold(5)]
    assert len(kfold) == 5
    i, j = kfold[0]
    assert (i != j).all()
