#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import warnings
import tempfile
import os

import pandas as pd

from fuel.datasets import H5PYDataset
from fuel.utils import find_in_data_path
from fuel import config

class Dataset(H5PYDataset):

    """ Abstract base class providing an interface to the skchem data format."""

    def __init__(self, **kwargs):
        kwargs.setdefault('load_in_memory', True)
        super(Dataset, self).__init__(
            file_or_path=find_in_data_path(self.filename), **kwargs)

    @classmethod
    def load_set(cls, set_name, sources=()):

        """ Load the sources for a single set.

        Args:
            set_name (str):
                The set name.
            sources (tuple[str]):
                The sources to return data for.

        Returns:
            tuple[np.array]
                The requested sources for the requested set.
        """
        if set_name == 'all':
            set_name = cls.set_names
        else:
            set_name = (set_name,)
        if sources == 'all':
            sources = cls.sources_names
        return cls(which_sets=set_name, sources=sources, load_in_memory=True).data_sources

    @classmethod
    def load_data(cls, sets=(), sources=()):

        """ Load a set of sources.

        Args:
            sets (tuple[str]):
                The sets to return data for.
            sources:
                The sources to return data for.

        Example:
            (X_train, y_train), (X_test, y_test) = Dataset.load_data(sets=('train', 'test'), sources=('X', 'y'))
        """

        for set_name in sets:
            yield cls.load_set(set_name, sources)

    @classmethod
    def read_frame(cls, key, *args, **kwargs):

        """ Load a set of features from the dataset as a pandas object.

        Args:
            key (str):
                The HDF5 key for required data.  Typically, this will be one of

                - structure: for the raw molecules
                - smiles: for the smiles
                - features/{feat_name}: for the features
                - targets/{targ_name}: for the targets

        Returns:
            pd.Series or pd.DataFrame or pd.Panel
                The data as a dataframe.
        """

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            data = pd.read_hdf(find_in_data_path(cls.filename), key, *args, **kwargs)
        if isinstance(data, pd.Panel):
            data = data.transpose(2, 1, 0)
        return data

    @classmethod
    def download(cls, output_directory=None, download_directory=None):

        """ Download the dataset and convert it.

        Args:
            output_directory (str):
                The directory to save the data to. Defaults to the first
                directory in the fuel data path.

            download_directory (str):
                The directory to save the raw files to. Defaults to a temporary
                directory.

        Returns:
            str:
                The path of the downloaded and processed dataset.
        """

        if not output_directory:
            output_directory = config.config['data_path']['yaml'].split(':')[0]

        output_directory = os.path.expanduser(output_directory)

        if not download_directory:
            download_directory = tempfile.mkdtemp()

        cls.downloader.download(directory=download_directory)
        return cls.converter.convert(directory=download_directory,
                                     output_directory=output_directory)
