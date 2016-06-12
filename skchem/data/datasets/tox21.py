#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..converters.tox21 import Tox21Converter
from ..downloaders.tox21 import Tox21Downloader

class Tox21(Dataset):
    filename = 'tox21.h5'
    downloader = Tox21Downloader
    converter = Tox21Converter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    Tox21.download()
