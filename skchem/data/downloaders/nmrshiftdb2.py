#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class NMRShiftDB2Downloader(Downloader):
    urls = ['https://sourceforge.net/p/nmrshiftdb2/code/HEAD/tree/trunk/snapshots/nmrshiftdb2withsignals.sd?format=raw']
    filenames = ['nmrshiftdb2.sdf']

if __name__ == '__main__':
    NMRShiftDB2Downloader.download(os.getcwd())
