#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import zipfile
import logging

LOGGER = logging.getLogger(__name__)

import numpy as np

from ... import io

from .base import Converter, default_pipeline, contiguous_order
from ...cross_validation import SimThresholdSplit

class BursiAmesConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='bursi_ames.h5'):

        zip_path = os.path.join(directory, 'cas_4337.zip')
        output_path = os.path.join(output_directory, output_filename)

        with zipfile.ZipFile(zip_path) as f:
            sdf_path = f.extract('cas_4337.sdf')

        data = io.read_sdf(sdf_path)
        data.index.name = 'batch'
        data['is_mutagen'] = (data['Ames test categorisation'] == 'mutagen').astype(np.uint8)
        ms, y = data.structure, data.is_mutagen
        pipeline = default_pipeline()
        ms, y = pipeline.transform_filter(ms, y)

        cv = SimThresholdSplit(ms,  min_threshold=0.6, n_jobs=-1)
        train, valid, test = cv.split((70, 15, 15))
        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))
        self.run(ms, y, output_path, splits=splits)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    LOGGER.info('Converting Bursi Ames Dataset...')
    BursiAmesConverter.convert()
