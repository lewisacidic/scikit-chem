#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.datasets.chembl

The ChEMBL dataset.
"""

import logging

from .base import Dataset
from ..converters.chembl import ChEMBLConverter
from ..downloaders.chembl import ChEMBLDownloader


class ChEMBL(Dataset):
    filename = 'chembl.h5'
    downloader = ChEMBLDownloader
    converter = ChEMBLConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    ChEMBL.download()
