#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.coverters.example

Formatter for the example dataset.
"""

import os

import pandas as pd
import numpy as np

from .base import Converter, contiguous_order, Feature
from ...pipeline import Pipeline
from ...io import read_sdf
from ...cross_validation import SimThresholdSplit
from ...descriptors import MorganFeaturizer
from ...standardizers import ChemAxonStandardizer

class DiversityConverter(Converter):

    """ Example Converter, using the NCI DTP Diversity Set III.  """

    def __init__(self, directory, output_directory, output_filename='diversity.h5'):

        output_path = os.path.join(output_directory, output_filename)

        ms = self.parse_file(os.path.join(directory, 'structures.sdf'))
        y = self.synthetic_targets(ms.index)

        pipeline = Pipeline([ChemAxonStandardizer(keep_failed=True)])

        cv = SimThresholdSplit(ms, min_threshold=0.6, n_jobs=-1)
        train, valid, test = cv.split((70, 15, 15))
        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))

        features = [Feature(fper=MorganFeaturizer(), key='X_morg', axis_names=['batch', 'features'])]

        self.run(ms, y, output_path, splits=splits, features=features)

    def parse_file(self, path):
        return read_sdf(path).structure

    def synthetic_targets(self, index):
        return pd.Series(np.random.choice([0, 1], size=len(index)), index=index)