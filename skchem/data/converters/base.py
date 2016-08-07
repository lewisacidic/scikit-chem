#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.converters.base

Defines the base converter class.
"""

import warnings
import logging
import os
from collections import namedtuple

import numpy as np
import pandas as pd
import h5py
from fuel.datasets import H5PYDataset

from ... import forcefields
from ... import filters
from ... import descriptors
from ... import standardizers
from ... import pipeline

logger = logging.getLogger(__name__)


def default_pipeline():
    """ Return a default pipeline to be used for general datasets. """
    return pipeline.Pipeline([
        standardizers.ChemAxonStandardizer(keep_failed=True, warn_on_fail=False),
        forcefields.UFF(add_hs=True, warn_on_fail=False),
        filters.OrganicFilter(),
        filters.AtomNumberFilter(above=5, below=100, include_hydrogens=True),
        filters.MassFilter(below=1000)
    ])

DEFAULT_PYTABLES_KW = {
    'complib': 'bzip2',
    'complevel': 9
}

def contiguous_order(to_order, splits):
    """ Determine a contiguous order from non-overlapping splits, and put data in that order.

    Args:
        to_order (iterable<pd.Series, pd.DataFrame, pd.Panel>):
            The pandas objects to put in contiguous order.
        splits (iterable<pd.Series>):
            The non-overlapping splits, as boolean masks.

    Returns:
        iterable<pd.Series, pd.DataFrame, pd.Panel>: The data in contiguous order.
    """

    member = pd.Series(0, index=splits[0].index)
    for i, split in enumerate(splits):
        member[split] = i
    idx = member.sort_values().index
    return (order.reindex(idx) for order in to_order)

Feature = namedtuple('Feature', ['fper', 'key', 'axis_names'])


def default_features():
    return (
        Feature(fper=descriptors.MorganFeaturizer(),
                key='X_morg',
                axis_names=['batch', 'features']),
        Feature(fper=descriptors.PhysicochemicalFeaturizer(),
                key='X_pc',
                axis_names=['batch', 'features']),
        Feature(fper=descriptors.AtomFeaturizer(max_atoms=100),
                key='A',
                axis_names=['batch', 'atom_idx', 'features']),
        Feature(fper=descriptors.GraphDistanceTransformer(max_atoms=100),
                key='G',
                axis_names=['batch', 'atom_idx', 'atom_idx']),
        Feature(fper=descriptors.SpacialDistanceTransformer(max_atoms=100),
                key='G_d',
                axis_names=['batch', 'atom_idx', 'atom_idx']),
        Feature(fper=descriptors.ChemAxonFeaturizer(features='all'),
                key='X_cx',
                axis_names=['batch', 'features']),
        Feature(fper=descriptors.ChemAxonAtomFeaturizer(features='all', max_atoms=100),
                key='A_cx',
                axis_names=['batch', 'atom_idx', 'features'])
    )


class Split(object):

    def __init__(self, mask, name, converter):
        self.mask = mask
        self.name = name
        self.converter = converter

    @property
    def contiguous(self):
        diff = np.ediff1d(self.mask.astype(int))
        if self.mask.iloc[0] != 0:
            diff[0] = 1
        if self.mask.iloc[-1] != 0:
            diff[-1] = -1
        return sum(diff == -1) == 1 or sum(diff == 1) == 1

    @property
    def indices(self):
        return np.nonzero(self.mask)[0]

    def save(self):
        self.converter.data_file[self.name + '_indices'] = self.indices
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.mask.to_hdf(self.converter.data_file.filename, '/indices/' + self.name)

    @property
    def ref(self):
        return self.converter.data_file[self.name + '_indices'].ref

    def to_dict(self):
        idx = self.indices
        if self.contiguous:
            low, high = min(idx), max(idx)
            return {source: (low, high) for source in self.converter.source_names}
        else:
            return {source: (-1, -1, self.ref) for source in self.converter.source_names}


class Converter(object):
    """ Create a fuel dataset from molecules and targets. """

    def __init__(self, directory, output_directory, output_filename='default.h5'):
        raise NotImplemented

    def run(self, ms, y, output_path, splits=None, features=None, pytables_kws=DEFAULT_PYTABLES_KW):

        """
           Args:
        ms (pd.Series):
            The molecules of the dataset.
        ys (pd.Series or pd.DataFrame):
            The target labels of the dataset.
        output_path (str):
            The path to which the dataset should be saved.
        features (list[Feature]):
            The features to calculate. Defaults are used if `None`.
        splits (iterable<(name, split)>):
            An iterable of name, split tuples. Splits are provided as boolean arrays of the whole data.
        """

        self.output_path = output_path
        self.pytables_kws = pytables_kws
        self.features = features if features is not None else default_features()
        self.feature_names = [feat.key for feat in self.features]
        self.task_names = ['y']
        self.splits = [Split(split, name, self) for name, split in splits]

        self.create_file(output_path)

        self.save_splits()
        self.save_molecules(ms)
        self.save_targets(y)
        self.save_features(ms)

    @property
    def source_names(self):
        return self.feature_names + self.task_names

    @property
    def split_names(self):
        return self.splits

    def create_file(self, path):
        logger.info('Creating h5 file at %s...', self.output_path)
        self.data_file = h5py.File(path, 'w')
        return self.data_file

    def save_molecules(self, mols):

        """ Save the molecules to the data file. """

        logger.info('Writing molecules to file...')
        logger.debug('Writing %s molecules to %s', len(mols), self.data_file.filename)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mols.to_hdf(self.data_file.filename, 'structure', **self.pytables_kws)
            mols.apply(lambda m: m.to_smiles().encode('utf-8')).to_hdf(self.data_file.filename, 'smiles')

    def save_frame(self, data, name, prefix='targets'):

        """ Save the a frame to the data file. """

        logger.info('Writing %s', name)
        logger.debug('Writing data of shape %s to %s', data.shape, self.data_file.filename)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if len(data.shape) > 2:
                data = data.transpose(2, 1, 0)  # panel serializes backwards for some reason...
            data.to_hdf(self.data_file.filename,
                        key='/{prefix}/{name}'.format(prefix=prefix, name=name),
                        **self.pytables_kws)

        if isinstance(data, pd.Series):
            self.data_file[name] = h5py.SoftLink('/{prefix}/{name}/values'.format(prefix=prefix, name=name))
            self.data_file[name].dims[0].label = data.index.name

        elif isinstance(data, pd.DataFrame):
            self.data_file[name] = h5py.SoftLink('/{prefix}/{name}/block0_values'.format(prefix=prefix, name=name))
            self.data_file[name].dims[0].label = data.index.name
            self.data_file[name].dims[1].label = data.columns.name

        elif isinstance(data, pd.Panel):
            self.data_file[name] = h5py.SoftLink('/{prefix}/{name}/block0_values'.format(prefix=prefix, name=name))
            self.data_file[name].dims[0].label = data.minor_axis.name # as panel serializes backwards
            self.data_file[name].dims[1].label = data.major_axis.name
            self.data_file[name].dims[2].label = data.items.name

    def save_targets(self, y):

        self.save_frame(y, name='y', prefix='targets')

    def save_features(self, ms):

        """ Save all features for the dataset. """
        logger.debug('Saving features')
        for feat in self.features:
            self._save_feature(ms, feat)

    def _save_feature(self, ms, feat):

        """ Calculate and save a feature to the data file. """
        logger.info('Calculating %s', feat.key)

        fps = feat.fper.transform(ms)
        self.save_frame(fps, name=feat.key, prefix='feats')

    def save_splits(self):

        """ Save the splits to the data file. """

        logger.info('Producing dataset splits...')
        for split in self.splits:
            split.save()
        split_dict = {split.name: split.to_dict() for split in self.splits}
        splits = H5PYDataset.create_split_array(split_dict)
        logger.debug('split: %s', splits)
        logger.info('Saving splits...')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.data_file.attrs['split'] = splits

    @classmethod
    def convert(cls, **kwargs):
        kwargs.setdefault('directory', os.getcwd())
        kwargs.setdefault('output_directory', os.getcwd())

        return cls(**kwargs).output_path,

    @classmethod
    def fill_subparser(cls, subparser):
        return cls.convert
