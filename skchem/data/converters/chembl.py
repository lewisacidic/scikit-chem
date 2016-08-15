#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.converters.chembl

Dataset constructor for ChEMBL
"""
import logging
import pandas as pd
import os

from .base import Converter, default_pipeline, contiguous_order, Feature
from ...cross_validation import SimThresholdSplit
from ... import descriptors

LOGGER = logging.getLogger(__name__)


class ChEMBLConverter(Converter):

    """ Converter for the ChEMBL dataset. """

    def __init__(self, directory, output_directory, output_filename='chembl.h5'):

        output_path = os.path.join(output_directory, output_filename)

        infile = os.path.join(directory, 'chembl_raw.h5')
        ms, y = self.parse_infile(infile)

        pipeline = default_pipeline()

        ms, y = pipeline.transform_filter(ms, y)

        cv = SimThresholdSplit(min_threshold=0.6, n_jobs=-1).fit(ms)
        train, valid, test = cv.split((70, 15, 15))
        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))

        features = (
        Feature(fper=descriptors.MorganFeaturizer(),
                key='X_morg',
                axis_names=['batch', 'features']),
        Feature(fper=descriptors.PhysicochemicalFeaturizer(),
                key='X_pc',
                axis_names=['batch', 'features']),
        Feature(fper=descriptors.AtomFeaturizer(max_atoms=100),
                key='A',
                axis_names=['batch', 'atom_idx', 'features']),
        Feature(fper=descriptors.GraphDistanceTransformer(max_atoms=100),
                key='G',
                axis_names=['batch', 'atom_idx', 'atom_idx']),
        Feature(fper=descriptors.SpacialDistanceTransformer(max_atoms=100),
                key='G_d'))

        self.run(ms, y, output_path, splits=splits)


    def parse_infile(self, filename):

        ms = pd.read_hdf(filename, 'structure')
        y = pd.read_hdf(filename, 'targets/Y')
        return ms, y

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    LOGGER.info('Converting ChEMBL...')
    ChEMBLConverter.convert()
