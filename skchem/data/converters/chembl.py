#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.converters.chembl

Dataset constructor for ChEMBL
"""

import pandas as pd
import os

from .base import Converter, default_pipeline, contiguous_order
from ...cross_validation import SimThresholdSplit

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
        self.run(ms, y, output_path, splits=splits)


    def parse_infile(self, filename):

        ms = pd.read_hdf(filename, 'structure')
        y = pd.read_hdf(filename, 'targets/Y')
        return ms, y