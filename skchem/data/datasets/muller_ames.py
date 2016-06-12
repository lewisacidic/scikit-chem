#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..converters.muller_ames import MullerAmesConverter
from ..downloaders.muller_ames import MullerAmesDownloader

class MullerAmes(Dataset):
    filename = 'muller_ames.h5'
    downloader = MullerAmesDownloader
    converter = MullerAmesConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    MullerAmes.download()
