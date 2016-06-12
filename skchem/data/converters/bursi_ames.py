#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import zipfile
import logging

logger = logging.getLogger(__name__)

import numpy as np

from ... import io

from .base import Converter

class BursiAmesConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='bursi_ames.h5'):

        zip_path = os.path.join(directory, 'cas_4337.zip')
        output_path = os.path.join(output_directory, output_filename)

        with zipfile.ZipFile(zip_path) as f:
            sdf_path = f.extract('cas_4337.sdf')

        data = io.read_sdf(sdf_path)
        data['is_mutagen'] = (data['Ames test categorisation'] == 'mutagen').astype(np.uint8)

        data = self.standardize(data)
        data = self.filter(data)

        ms, y = data.structure, data.is_mutagen
        self.run(ms, y, output_path, contiguous=True)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    BursiAmesConverter.convert()
