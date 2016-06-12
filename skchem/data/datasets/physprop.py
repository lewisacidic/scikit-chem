#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..converters.physprop import PhysPropConverter
from ..downloaders.physprop import PhysPropDownloader

class PhysProp(Dataset):
    filename = 'physprop.h5'
    downloader = PhysPropDownloader
    converter = PhysPropConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    PhysProp.download()
