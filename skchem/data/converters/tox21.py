#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.data.transformers.tox21

Module defining transformation techniques for tox21.
"""

import zipfile
import os
import logging
LOGGER = logging.getLogger(__name__)

import numpy as np
import pandas as pd

from .base import Converter, default_pipeline
from ... import io
from ... import core

class Tox21Converter(Converter):

    """ Class to build tox21 dataset.

    """
    def __init__(self, directory, output_directory, output_filename='tox21.h5'):

        output_path = os.path.join(output_directory, output_filename)

        # extract data
        train, valid, test = self.extract(directory)

        # read data
        train = self.read_train(train)
        valid = self.read_valid(valid)
        test = self.read_test(test, os.path.join(directory, 'test.txt'))

        # combine into full dataset
        data = pd.concat([train, valid, test], keys=['train', 'valid', 'test']).sort_index()
        data.index.names = 'ds', 'id'

        ms, y = data.structure, data.drop('structure', axis=1)

        pipeline = default_pipeline()
        ms, y = pipeline.transform_filter(ms, y)

        # generate splits
        ms, y = ms.reset_index(0), y.reset_index(0)
        split_arr = ms.pop('ds')
        y.pop('ds')

        splits = [(split, split_arr == split) for split in ('train', 'valid', 'test')]

        y.columns.name = 'tasks'

        # call the Converter to make the final dataset
        self.run(ms, y, output_path, splits=splits)

    @staticmethod
    def fix_id(s):
        return s.split('-')[0]

    @staticmethod
    def fix_assay_name(s):
        return s.replace('-', '_')

    @staticmethod
    def patch_test(test):
        test_1 = pd.Series({
            'structure': core.Mol.from_smiles('FC(F)(F)c1[nH]c(c(C#N)c1Br)C1=CC=C(Cl)C=C1', name='NCGC00357062'),
            'stochiometry': 0,
            'Compound ID': 'NCGC00357062',
            'Sample ID': 'NCGC00357062-01'}, name='NCGC00357062')
        test['NCGC00357062'] = test_1
        return test

    def read_train(self, train):

        train = io.read_sdf(train)
        train.columns = train.columns.to_series().apply(self.fix_assay_name)
        train.index = train.index.to_series().apply(self.fix_id)
        self.assays = train.columns[-12:]
        self.keep_cols = ['structure'] + self.assays.tolist()
        train[self.assays] = train[self.assays].astype(float)
        train = train[self.keep_cols]
        train = train.sort_index()
        ms = train.structure[~train.index.duplicated()]
        train = train[self.assays].groupby(train.index).max()
        train = ms.to_frame().join(train)
        return train

    def read_valid(self, valid):

        valid = io.read_sdf(valid)
        valid.columns = valid.columns.to_series().apply(self.fix_assay_name)
        valid = valid[self.keep_cols]
        valid[self.assays] = valid[self.assays].astype(float)
        return valid

    def read_test(self, test, test_data):

        test = io.read_sdf(test)
        test = self.patch_test(test)
        test_data = pd.read_table(test_data)
        test_data['Sample ID'] = test_data['Sample ID'].apply(self.fix_id)
        test = test.join(test_data.set_index('Sample ID'))

        test.columns = test.columns.to_series().apply(self.fix_assay_name)
        test = test[self.keep_cols]
        test[test == 'x'] = np.nan
        test[self.assays] = test[self.assays].astype(float)
        return test

    def extract(self, directory):

        with zipfile.ZipFile(os.path.join(directory, 'train.sdf.zip')) as f:
            train = f.extract('tox21_10k_data_all.sdf')

        with zipfile.ZipFile(os.path.join(directory, 'valid.sdf.zip')) as f:
            valid = f.extract('tox21_10k_challenge_test.sdf')

        with zipfile.ZipFile(os.path.join(directory, 'test.sdf.zip')) as f:
            test = f.extract('tox21_10k_challenge_score.sdf')

        return train, valid, test

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    LOGGER.info('Converting Tox21 Dataset...')
    Tox21Converter.convert()
