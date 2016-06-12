#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..converters.bursi_ames import BursiAmesConverter
from ..downloaders.bursi_ames import BursiAmesDownloader

class BursiAmes(Dataset):
    filename = 'bursi_ames.h5'
    downloader = BursiAmesDownloader
    converter = BursiAmesConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    BursiAmes.download()
