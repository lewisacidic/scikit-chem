#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class BursiAmesDownloader(Downloader):
    urls = ['http://cheminformatics.org/datasets/bursi/cas_4337.zip']
    filenames = ['cas_4337.zip']

if __name__ == '__main__':
    BursiAmesDownloader.download(os.getcwd())
