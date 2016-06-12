#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import logging

from .base import Dataset
from ..downloaders.bradley_open_mp import BradleyOpenMPDownloader
from ..converters.bradley_open_mp import BradleyOpenMPConverter

class BradleyOpenMP(Dataset):
    filename = 'bradley_open_mp.h5'
    downloader = BradleyOpenMPDownloader
    converter = BradleyOpenMPConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    BradleyOpenMP.download()
