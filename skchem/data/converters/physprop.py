#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
import zipfile
import logging
LOGGER = logging.getLogger(__name__)

import pandas as pd
import numpy as np

from ... import io
from .base import Converter, contiguous_order

from ...cross_validation import SimThresholdSplit

TXT_COLUMNS = [l.lower() for l in """CAS
Formula
Mol_Weight
Chemical_Name
WS
WS_temp
WS_type
WS_reference
LogP
LogP_temp
LogP_type
LogP_reference
VP
VP_temp
VP_type
VP_reference
DC_pKa
DC_temp
DC_type
DC_reference
henry_law Constant
HL_temp
HL_type
HL_reference
OH
OH_temp
OH_type
OH_reference
BP_pressure
MP
BP
FP""".split('\n')]

class PhysPropConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='physprop.h5'):

        output_path = os.path.join(output_directory, output_filename)

        sdf, txt = self.extract(directory)
        mols, data = self.process_sdf(sdf), self.process_txt(txt)

        LOGGER.debug('Compounds with data extracted: %s', len(data))

        data = mols.to_frame().join(data)
        data = self.drop_inconsistencies(data)

        y = self.process_targets(data)
        LOGGER.debug('Compounds with experimental: %s', len(y))

        data = data.ix[y.index]
        data.columns.name = 'targets'
        ms, y = data.structure, data.drop('structure', axis=1)

        cv = SimThresholdSplit(ms, min_threshold=0.6, block_width=4000, n_jobs=-1)
        train, valid, test = cv.split((70, 15, 15))

        (ms, y, train, valid, test) = contiguous_order((ms, y, train, valid, test), (train, valid, test))
        splits = (('train', train), ('valid', valid), ('test', test))

        self.run(ms, y, output_path=output_path, splits=splits)

    def extract(self, directory):
        LOGGER.info('Extracting from %s', directory)
        with zipfile.ZipFile(os.path.join(directory, 'phys_sdf.zip')) as f:
            sdf = f.extract('PhysProp.sdf')
        with zipfile.ZipFile(os.path.join(directory, 'phys_txt.zip')) as f:
            txt = f.extract('PhysProp.txt')
        return sdf, txt

    def process_sdf(self, path):
        LOGGER.info('Processing sdf at %s', path)
        mols = io.read_sdf(path, read_props=False).structure
        mols.index = mols.apply(lambda m: m.GetProp('CAS'))
        mols.index.name = 'cas'
        LOGGER.debug('Structures extracted: %s', len(mols))
        return mols

    def process_txt(self, path):
        LOGGER.info('Processing txt at %s', path)
        data = pd.read_table(path, header=None, engine='python').iloc[:, :32]
        data.columns = TXT_COLUMNS
        data_types = data.columns[[s.endswith('_type') for s in data.columns]]
        data[data_types] = data[data_types].fillna('NAN')
        data = data.set_index('cas')
        return data

    def drop_inconsistencies(self, data):
        LOGGER.info('Dropping inconsistent data...')
        formula = data.structure.apply(lambda m: m.to_formula())
        LOGGER.info('Inconsistent compounds: %s', (formula != data.formula).sum())
        data = data[formula == data.formula]
        return data

    def process_targets(self, data):
        LOGGER.info('Dropping estimated data...')
        data = pd.concat([self.process_logS(data),
                          self.process_logP(data),
                          self.process_mp(data),
                          self.process_bp(data)], axis=1)
        LOGGER.info('Dropped compounds: %s', data.isnull().all(axis=1).sum())
        data = data[data.notnull().any(axis=1)]
        LOGGER.debug('Compounds with experimental activities: %s', len(data))
        return data

    def process_logS(self, data):
        cleaned = pd.DataFrame(index=data.index)
        S = 0.001 * data.ws / data.mol_weight
        logS = np.log10(S)
        return logS[data.ws_type == 'EXP']

    def process_logP(self, data):
        logP = data.logp[data.logp_type == 'EXP']
        return logP[logP > -10]

    def process_mp(self, data):
        return data.mp.apply(self.fix_temp)

    def process_bp(self, data):
        return data.bp.apply(self.fix_temp)

    @staticmethod
    def fix_temp(s, mean_range=5):
        try:
            return float(s)
        except ValueError:
            if '<' in s or '>' in s:
                return np.nan
            s = s.strip(' dec')
            s = s.strip(' sub')
            if '-' in s and mean_range:
                rng = [float(n) for n in s.split('-')]
                if len(rng) > 2:
                    return np.nan
                if np.abs(rng[1] - rng[0]) < mean_range:
                    return (rng[0] + rng[1])/2
            try:
                return float(s)
            except ValueError:
                return np.nan



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    LOGGER.info('Converting PhysProp Dataset...')
    PhysPropConverter.convert()
