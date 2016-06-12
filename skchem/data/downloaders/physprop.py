#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class PhysPropDownloader(Downloader):
    filenames = ['phys_sdf.zip', 'phys_txt.zip']
    urls = ['http://esc.syrres.com/interkow/Download/' + f for f in filenames]

if __name__ == '__main__':
    PhysPropDownloader.download(os.getcwd())
