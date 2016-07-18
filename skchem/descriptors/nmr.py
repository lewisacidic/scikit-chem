#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.nmr

Module for NMR prediction.
"""
import tempfile
import subprocess
import re
import logging
import time

import numpy as np
import pandas as pd

from ..io import write_sdf
from ..utils import NamedProgressBar
from .. import core

LOGGER = logging.getLogger(__name__)

def mol_count(filename):
    return sum(1 for l in open(filename, 'rb') if l == b'##PEAKASSIGNMENTS=(XYMA)\r\n')


class ChemAxonNMRPredictor(object):

    def __init__(self, element='C', max_atoms='auto'):
        self.element = element
        self._detect_max_atoms = False
        self.max_atoms = max_atoms

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, val):
        val = val.upper()
        if val in ('C', 'H'):
            self._element = val
        else:
            raise ValueError('ChemAxon can only predict 1H or 13C shifts.')

    @property
    def max_atoms(self):
        return self._max_atoms

    @max_atoms.setter
    def max_atoms(self, arg):
        if arg == 'auto':
            self._detect_max_atoms = True
            self._max_atoms = -1
        else:
            self._max_atoms = arg

    @property
    def index(self):
        return pd.RangeIndex(self.max_atoms, name='atom_idx')

    @property
    def feature_names(self):
        return self.index

    def transform(self, obj):
        if isinstance(obj, core.Mol):
            return self._transform_mol(obj)
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_mol(self, mol):
        # make into series then use self._transform_mol
        ser = pd.Series([mol], name=mol.name)
        res = self._transform_series(ser)
        return res.iloc[0]

    def _transform_series(self, ser):

        if self._detect_max_atoms:
            self.max_atoms = ser.ms.atoms.apply(len).max()

        with tempfile.NamedTemporaryFile(suffix='.sdf') as infile, tempfile.NamedTemporaryFile() as outfile:

            write_sdf(ser, infile.name)
            args = ['cxcalc', infile.name, '-o', outfile.name] + [self.element.lower() + 'nmr']

            LOGGER.info('Running: ' + ' '.join(args))

            p = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            bar = NamedProgressBar(name=self.__class__.__name__, max_value=len(ser))

            while p.poll() is None:
                bar.update(mol_count(outfile.name))
                time.sleep(1)
            bar.update(len(ser))
            p.wait()

            res = self.parse_output(outfile.name, n_mols=len(ser))
            return pd.DataFrame(res, index=ser.index, columns=self.index)


    def parse_output(self, outfile, n_mols=None):

        """ Read the NMR output file. """

        if n_mols is None:
            n_mols = mol_count(outfile)

        res = np.repeat(np.nan, n_mols * self.max_atoms).reshape(n_mols, self.max_atoms)
        regex = re.compile(b'\((-?\d+.\d+),\d+,[A-Z],<([0-9\,]+)>\)\r\n')

        mol_idx = 0

        with open(outfile, 'rb') as f:
            # loop through the file - inner loop will also advance the pointer
            for l in f:
                if l == b'##PEAKASSIGNMENTS=(XYMA)\r\n':
                    for row in f:
                        if row == b'##END=\r\n':
                            break
                        else:
                            LOGGER.debug('Row to parse: %s', row)
                            shift, idxs = regex.match(row).groups()
                            shift, idxs = float(shift), [int(idx) for idx in idxs.split(b',')]
                            for atom_idx in idxs:
                                res[mol_idx, atom_idx] = shift
                    mol_idx += 1
        return res




