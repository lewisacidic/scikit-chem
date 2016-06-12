#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import warnings
import logging
import os
import functools
from collections import namedtuple

import numpy as np
import pandas as pd

import h5py
from fuel.datasets import H5PYDataset

from ... import filters
from ... import descriptors
from ... import cross_validation
from ... import standardizers

logger = logging.getLogger(__name__)


Feature = namedtuple('Feature', ['fper', 'key', 'axis_names'])

DEFAULT_FEATURES = (
    Feature(fper=descriptors.MorganFingerprinter(),
            key='X_morg',
            axis_names=['batch', 'features']),
    Feature(fper=descriptors.PhysicochemicalFingerprinter(),
            key='X_pc',
            axis_names=['batch', 'features']),
    Feature(fper=descriptors.AtomFeatureCalculator(),
            key='A',
            axis_names=['batch', 'atom_idx', 'features']),
    Feature(fper=descriptors.GraphDistanceCalculator(),
            key='G',
            axis_names=['batch', 'atom_idx', 'atom_idx']))

Filter = namedtuple('Filter', ['filter', 'kwargs'])

DEFAULT_FILTERS = (
    Filter(filters.is_organic, {}),
    Filter(filters.n_atoms, {'above': 5, 'below': 75}),
    Filter(filters.mass, {'below': 1000})
)

DEFAULT_STANDARDIZER = standardizers.ChemAxonStandardizer(keep_failed=True)

class Converter(object):
    """ Create a fuel dataset from molecules and targets.

    Args:
        ms (pd.Series):
            The molecules of the dataset.
        ys (pd.Series or pd.DataFrame):
            The target labels of the dataset.
        output_path (str):
            The path to which the dataset should be saved.
        features (list[Feature]):
            The features to calculate. Defaults are provided.
        splits (dict):
            A dictionary of different splits provided.
            The keys should be the split name, and values an array of indices.
            Alternatively, if `contiguous_splits` is `True`, the keys should be
            the split name, and the values a tuple of start and stop.
            If `None`, use `skchem.cross_validation.SimThresholdSplit`
    """


    def __init__(self, directory, output_directory, output_filename='default.h5'):
        raise NotImplemented

    def run(self, ms, y, output_path,
                features=DEFAULT_FEATURES, splits=None, contiguous=False):

        self.contiguous = contiguous
        self.output_path = output_path
        self.features = features
        self.feature_names = [feat.key for feat in self.features] + ['y']

        self.create_file(output_path)

        if not splits:
            splits, idx = self.create_splits(ms)
            ms, y = ms.ix[idx], y.ix[idx]

        split_dict = self.process_splits(splits)

        self.save_splits(split_dict)
        self.save_molecules(ms)
        self.save_targets(y)
        self.save_features(ms)

    def create_file(self, path):
        logger.info('Creating h5 file at %s...', self.output_path)
        self.data_file = h5py.File(path, 'w')
        return self.data_file

    def filter(self, data, filters=DEFAULT_FILTERS):

        """ Filter the compounds according to the usual filters. """
        logger.info('Filtering %s compounds', len(data))
        if isinstance(data, pd.DataFrame):
            ms = data.structure
        else:
            ms = data
        filt = functools.reduce(lambda a, b: a & b, (ms.apply(filt.filter, **filt.kwargs) for filt in filters))
        logger.info('Filtered out %s compounds', (~filt).sum())

        return data[filt]


    def standardize(self, data, standardizer=DEFAULT_STANDARDIZER):

        """ Standardize the compounds. """
        logger.info('Standardizing %s compounds', len(data))
        return standardizer.transform(data)


    def save_molecules(self, mols):

        """ Save the molecules to the data file. """

        logger.info('Writing molecules to file...')
        logger.debug('Writing %s molecules to %s', len(mols), self.data_file.filename)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mols.to_hdf(self.data_file.filename, 'structure')
            mols.apply(lambda m: m.to_smiles().encode('utf-8')).to_hdf(self.data_file.filename, 'smiles')

    def save_targets(self, y):

        """ Save the targets to the data file. """
        y_name = getattr(y, 'name', None)
        if not y_name:
            y_name = getattr(y.columns, 'name', None)
        if not y_name:
            y_name = 'targets'

        logger.info('Writing %s', y_name)
        logger.debug('Writing targets of shape %s to %s', y.shape, self.data_file.filename)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            y.to_hdf(self.data_file.filename, '/targets/' + y_name)

        if isinstance(y, pd.Series):
            self.data_file['y'] = h5py.SoftLink('/targets/{}/values'.format(y_name))
            self.data_file['y'].dims[0].label = 'batch'

        elif isinstance(y, pd.DataFrame):
            self.data_file['y'] = h5py.SoftLink('/targets/{}/block0_values'.format(y_name))
            self.data_file['y'].dims[0].label = 'batch'
            self.data_file['y'].dims[0].label = 'task'

    def save_features(self, ms):

        """ Save all features for the dataset. """
        logger.debug('Saving features')
        for feat in self.features:
            self._save_feature(ms, feat)

    def _save_feature(self, ms, feat):

        """ Calculate and save a feature to the data file. """
        logger.info('Calculating %s', feat.key)

        fps = feat.fper.transform(ms)
        if len(feat.axis_names) > 2:
            fps = fps.transpose(2, 1, 0) # panel serialize backwards for some reason...
        logger.debug('Writing features with shape %s to %s', fps.shape, self.data_file.filename)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fps.to_hdf(self.data_file.filename, 'features/{}'.format(feat.key))
        self.data_file[feat.key] = h5py.SoftLink('/features/{}/block0_values'.format(feat.key))
        self.data_file[feat.key].dims[0].label = feat.axis_names[0]
        self.data_file[feat.key].dims[1].label = feat.axis_names[1]
        if len(feat.axis_names) > 2:
            self.data_file[feat.key].dims[2].label = feat.axis_names[2]

    def create_splits(self, ms, contiguous=True):

        """ Create a split dict for fuel from mols, using SimThresholdSplit.

        Args:
            ms (pd.Series):
                The molecules to use to design the splits.
            contiguous (bool):
                Whether the split should be contiguous.  This allows for more
                efficient loading times.  This usually is the appropriate if
                there are no other splits for the dataset, and will reorder
                the dataset.
        Returns:
            (dict, idx)
                The split dict, and the index to align the data with.
        """

        logger.info('Creating Similarity Threshold splits...')
        cv = cross_validation.SimThresholdSplit(ms, memory_optimized=True)
        train, valid, test = cv.split((70, 15, 15))

        def bool_to_index(ser):
            return np.nonzero(ser.values)[0]

        if self.contiguous:
            dset = pd.Series(0, ms.index)
            dset[train] = 0
            dset[valid] = 1
            dset[test] = 2
            dset = dset.sort_values()
            idx = dset.index
            train_split = bool_to_index(dset == 0)
            valid_split = bool_to_index(dset == 1)
            test_split = bool_to_index(dset == 2)
            print('train', train_split)
            print('valid', valid_split)
            print('test', test_split)
            def min_max(split):
                return min(split), max(split)

            splits = {
                'train': min_max(train_split),
                'valid': min_max(valid_split),
                'test': min_max(test_split)
            }

        else:

            idx = ms.index

            splits = {
                'train': bool_to_index(train),
                'valid': bool_to_index(valid),
                'test': bool_to_index(test)
            }

        return splits, idx

    def process_splits(self, splits, contiguous=False):

        """ Create a split dict for fuel from provided indexes. """

        logger.info('Creating split array.')

        split_dict = {}

        if self.contiguous:
            logger.debug('Contiguous splits.')
            for split_name, (start, stop) in splits.items():
                split_dict[split_name] = {feat: (start, stop, h5py.Reference()) for feat in self.feature_names}
        else:
            for split_name, split in splits.items():
                split_indices_name = '{}_indices'.format(split_name).encode('utf-8')
                logger.debug('Saving %s to %s', split_indices_name, self.data_file.filename)
                self.data_file[split_indices_name] = split
                split_ref = self.data_file[split_indices_name].ref
                split_dict[split_name] = {feat: (-1, -1, split_ref) for feat in self.feature_names}

        return split_dict

    def save_splits(self, split_dict):

        """ Save the splits to the data file. """

        logger.info('Producing dataset splits...')
        split = H5PYDataset.create_split_array(split_dict)
        logger.debug('split: %s', split)
        logger.info('Saving splits...')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.data_file.attrs['split'] = split

    @classmethod
    def convert(cls, **kwargs):
        kwargs.setdefault('directory', os.getcwd())
        kwargs.setdefault('output_directory', os.getcwd())

        return cls(**kwargs).output_path,

    @classmethod
    def fill_subparser(cls, subparser):
        return cls.convert
