#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..converters.nmrshiftdb2 import NMRShiftDB2Converter
from ..downloaders.nmrshiftdb2 import NMRShiftDB2Downloader

class NMRShiftDB2(Dataset):
    filename = 'nmrshiftdb2.h5'
    downloader = NMRShiftDB2Downloader
    converter = NMRShiftDB2Converter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    NMRShiftDB2.download()
