#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import logging
logger = logging.getLogger(__name__)

import pandas as pd

from ... import core
from .base import Converter

class BradleyOpenMPConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='bradley_open_mp.h5'):

        output_path = os.path.join(output_directory, output_filename)
        data = self.parse_data(os.path.join(directory, 'bradley_melting_point_dataset.xlsx'))
        data = self.filter_bad(data)
        data['structure'] = self.standardize(data.smiles)
        data = self.filter(data)
        ms, y = data.structure, self.fix_mp(data)
        self.run(ms, y, output_path=output_path)

    def parse_data(self, path):
        logger.info('Parsing data at %s...', path)
        return pd.read_excel(path, index_col=0)

    def filter_bad(self, data):
        logger.info('Removing manually annotated errors...')
        bad_data = data.donotuse.notnull()
        logger.debug('Removed %s', bad_data.sum())
        return data[~bad_data]

    def fix_mp(self, data):
        logger.info('Converting temperature to Kelvin...')
        return data.mpC + 278.15

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    LOGGER.info('Converting Bradley Open Melting Point Dataset...')
    BradleyOpenMPConverter.convert()
