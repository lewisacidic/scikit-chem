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
import skchem

from .base import Converter

from ... import standardizers

PATCHES = {
    '820-75-7': r'NNC(=O)CNC(=O)C=[N+]=[N-]',
    '2435-76-9': r'[N-]=[N+]=C1C=NC(=O)NC1=O',
    '817-99-2': r'NC(=O)CNC(=O)\C=[N+]=[N-]',
    '116539-70-9': r'CCCCN(CC(O)C1=C\C(=[N+]=[N-])\C(=O)C=C1)N=O',
    '115-02-6': r'NC(COC(=O)\C=[N+]=[N-])C(=O)O',
    '122341-55-3': r'NC(COC(=O)\C=[N+]=[N-])C(=O)O'
}

class MullerAmesConverter(Converter):

    def __init__(self, directory, output_directory, output_filename='muller_ames.h5'):

        """
        Args:
            directory (str):
                Directory in which input files reside.
            output_directory (str):
                Directory in which to save the converted dataset.
            output_filename (str):
                Name of the saved dataset. Defaults to `muller_ames.h5`.

        Returns:
            tuple of str:
                Single-element tuple containing the path to the converted dataset.
        """

        zip_path = os.path.join(directory, 'ci900161g_si_001.zip')
        output_path = os.path.join(output_directory, output_filename)

        with zipfile.ZipFile(zip_path) as f:
            f.extractall()

        # create dataframe
        data = pd.read_csv(os.path.join(directory, 'smiles_cas_N6512.smi'),
                           delimiter='\t', index_col=1,
                           converters={1: lambda s: s.strip()},
                           header=None, names=['structure', 'id', 'is_mutagen'])

        data = self.patch_data(data, PATCHES)

        data['structure'] = data.structure.apply(skchem.Mol.from_smiles)

        data = self.standardize(data)
        data = self.optimize(data)
        keep = self.filter(data)

        ms, ys = keep.structure, keep.is_mutagen

        indices = data.reset_index().index.difference(keep.reset_index().index)

        train = self.parse_splits(os.path.join('splits_train_N6512.csv'))
        train = self.drop_indices(train, indices)
        splits = self.create_split_dict(train, 'train')

        test = self.parse_splits(os.path.join(directory, 'splits_test_N6512.csv'))
        test = self.drop_indices(test, indices)
        splits.update(self.create_split_dict(test, 'test'))

        self.run(ms, ys, output_path, splits=splits)

    def patch_data(self, data, patches):
        """ Patch smiles in a DataFrame with rewritten ones that specify diazo
        groups in rdkit friendly way. """

        LOGGER.info('Patching data...')
        for cas, smiles in patches.items():
            data.loc[cas, 'structure'] = smiles

        return data

    def parse_splits(self, f_path):
        LOGGER.info('Parsing splits...')
        with open(f_path) as f:
            splits = [split for split in f.read().strip().splitlines()]

        splits = [[n for n in split.strip().split(',')] for split in splits]
        splits = [sorted(int(n) for n in split) for split in splits] # sorted ints
        return [np.array(split) - 1 for split in splits] # zero based indexing

    def drop_indices(self, splits, indices):
        LOGGER.info('Dropping failed compounds from split indices...')
        for i, split in enumerate(splits):
            split = split - sum(split > ix for ix in indices)
            splits[i] = np.delete(split, indices)

        return splits

    def create_split_dict(self, splits, name):
        return {'{}_{}'.format(name, i + 1): split \
                        for i, split in enumerate(splits)}

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    LOGGER.info('Converting Muller Ames Dataset...')
    MullerAmesConverter.convert()
