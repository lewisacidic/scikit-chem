#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import logging
logger = logging.getLogger(__name__)

import pandas as pd

from .base import Converter, default_pipeline, contiguous_order
from ...core import Mol
from ...cross_validation import SimThresholdSplit

class BradleyOpenMPConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='bradley_open_mp.h5'):

        output_path = os.path.join(output_directory, output_filename)
        data = self.parse_data(os.path.join(directory, 'bradley_melting_point_dataset.xlsx'))
        data = self.filter_bad(data)

        def parse_smiles(smi):
            try:
                return Mol.from_smiles(smi)
            except ValueError:
                return None

        data['structure'] = data.smiles.apply(parse_smiles)
        data = data[data.structure.notnull()]
        ms, y = data.structure, self.fix_mp(data)

        pipeline = default_pipeline()
        ms, y = pipeline.transform_filter(ms, y)

        cv = SimThresholdSplit(ms, min_threshold=0.6, n_jobs=-1)
        train, valid, test = cv.split((70, 15, 15))
        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))

        self.run(ms, y, output_path=output_path, splits=splits)

    @staticmethod
    def parse_data(path):
        logger.info('Parsing data at %s...', path)
        return pd.read_excel(path, index_col=0)

    @staticmethod
    def filter_bad(data):
        logger.info('Removing manually annotated errors...')
        bad_data = data.donotuse.notnull()
        logger.debug('Removed %s', bad_data.sum())
        return data[~bad_data]

    @staticmethod
    def fix_mp(data):
        logger.info('Converting temperature to Kelvin...')
        return data.mpC + 278.15

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    LOGGER.info('Converting Bradley Open Melting Point Dataset...')
    BradleyOpenMPConverter.convert()
