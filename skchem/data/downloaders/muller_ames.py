#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class MullerAmesDownloader(Downloader):
    urls = ['https://ndownloader.figshare.com/files/4523278']
    filenames = ['ci900161g_si_001.zip']

if __name__ == '__main__':
    MullerAmesDownloader.download(os.getcwd())
