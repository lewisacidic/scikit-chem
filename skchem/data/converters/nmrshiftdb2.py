#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import logging
import itertools
from collections import defaultdict

import pandas as pd
import numpy as np
from sklearn import metrics

from .base import Converter, default_pipeline, contiguous_order
from ... import io
from ... import utils
from ...cross_validation import SimThresholdSplit

LOGGER = logging.getLogger(__file__)

class NMRShiftDB2Converter(Converter):

    def __init__(self, directory, output_directory, output_filename='nmrshiftdb2.h5'):

        output_path = os.path.join(output_directory, output_filename)
        input_path = os.path.join(directory, 'nmrshiftdb2.sdf')
        data = self.parse_data(input_path)

        ys = self.get_spectra(data)
        ys = self.process_spectra(ys)
        ys = self.combine_duplicates(ys)
        self.log_dists(ys)
        self.log_duplicates(ys)
        ys = self.squash_duplicates(ys)

        c13s = self.to_frame(ys.loc[ys['13c'].notnull(), '13c'])
        data = data[['structure']].join(c13s, how='right')

        ms, y = data.structure, data.drop('structure', axis=1)
        pipeline = default_pipeline()
        ms, y = pipeline.transform_filter(ms, y)
        y.columns.name = 'shifts'

        cv = SimThresholdSplit(ms, min_threshold=0.6, block_width=4000, n_jobs=-1)
        train, valid, test = cv.split((70, 15, 15))

        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))

        self.run(ms, y, output_path=output_path, splits=splits)

    @staticmethod
    def parse_data(filepath):

        """ Reads the raw datafile. """

        LOGGER.info('Reading file: %s', filepath)
        data = io.read_sdf(filepath, removeHs=False, warn_bad_mol=False)
        data.index = data['nmrshiftdb2 ID'].astype(int)
        data.index.name = 'nmrshiftdb2_id'
        data.columns = data.columns.to_series().apply(utils.free_to_snail)
        data = data.sort_index()
        LOGGER.info('Read %s molecules.', len(data))
        return data

    @staticmethod
    def get_spectra(data):

        """ Retrieves spectra from raw data. """

        LOGGER.info('Retrieving spectra from raw data...')
        isotopes = [
            '1h',
            '11b',
            '13c',
            '15n',
            '17o',
            '19f',
            '29si',
            '31p',
            '33s',
            '73ge',
            '195pt'
        ]

        def is_spectrum(col_name, ele='c'):
            return any(isotope in col_name for isotope in isotopes)

        spectrum_cols = [c for c in data if is_spectrum(c)]
        data = data[spectrum_cols]

        def index_pair(s):
            return s[0], int(s[1])

        data.columns = pd.MultiIndex.from_tuples([index_pair(i.split('_')[1:]) for i in data.columns])
        return data

    @staticmethod
    def process_spectra(data):

        """ Turn the string representations found in sdf file into a dictionary. """

        def spectrum_dict(spectrum_string):
            if not isinstance(spectrum_string, str):
                return np.nan # no spectra are still nan
            if spectrum_string == '':
                return np.nan # empty spectra are nan
            sigs = spectrum_string.strip().strip('|').strip().split('|') # extract signals
            sig_tup = [tuple(s.split(';')) for s in sigs] # take tuples as (signal, coupling, atom)
            return {int(s[2]): float(s[0]) for s in sig_tup} # make spectrum a dictionary of atom to signal

        return data.applymap(spectrum_dict)

    @staticmethod
    def combine_duplicates(data):

        """ Collect duplicate spectra into one dictionary. All shifts are collected into lists. """

        def aggregate_dicts(ds):
            res = defaultdict(list)
            for d in ds:
                if not isinstance(d, dict): continue
                for k, v in d.items():
                    res[k].append(v)
            return dict(res) if len(res) else np.nan

        return data.groupby(level=0, axis=1).apply(lambda s: s.apply(aggregate_dicts, axis=1))

    @staticmethod
    def squash_duplicates(data):

        """ Take the mean of all the duplicates.  This is where we could do a bit more checking. """

        def squash(d):
            if not isinstance(d, dict):
                return np.nan
            else:
                return {k: np.mean(v) for k, v in d.items()}

        return data.applymap(squash)

    @staticmethod
    def to_frame(data):

        """ Convert a series of dictionaries to a dataframe. """
        res = pd.DataFrame(data.tolist(), index=data.index)
        res.columns.name = 'atom_idx'
        return res

    @staticmethod
    def extract_duplicates(data, kind='13c'):

        """ Get all 13c duplicates.  """

        def is_duplicate(ele):
            if not isinstance(ele, dict):
                return False
            else:
                return len(list(ele.values())[0]) > 1

        return data.loc[data[kind].apply(is_duplicate), kind]

    @staticmethod
    def log_dists(data):

        def n_spect(ele):
            return isinstance(ele, dict)

        def n_shifts(ele):
            return len(ele) if isinstance(ele, dict) else 0

        def log_message(func):
            return '  '.join('{k}: {v}'.format(k=k, v=v) for k, v in data.applymap(func).sum().to_dict().items())

        LOGGER.info('Number of spectra: %s', log_message(n_spect))
        LOGGER.info('Extracted shifts: %s', log_message(n_shifts))


    def log_duplicates(self, data):

        for kind in '1h', '13c':
            dups = self.extract_duplicates(data, kind)
            LOGGER.info('Number of duplicate %s spectra: %s', kind, len(dups))
            res = pd.DataFrame(sum((list(itertools.combinations(l, 2)) for s in dups for k, l in s.items()), []))
            LOGGER.info('Number of duplicate %s pairs: %f', kind, len(res))
            LOGGER.info('MAE for duplicate %s: %.4f', kind, metrics.mean_absolute_error(res[0], res[1]))
            LOGGER.info('MSE for duplicate %s: %.4f', kind, metrics.mean_squared_error(res[0], res[1]))
            LOGGER.info('r2 for duplicate %s: %.4f', kind, metrics.r2_score(res[0], res[1]))


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    LOGGER.info('Converting NMRShiftDB2 Dataset...')
    NMRShiftDB2Converter.convert()