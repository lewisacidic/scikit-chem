#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# file title

Description
"""

import os
from .base import Downloader

class DiversityDownloader(Downloader):
    urls = ['https://wiki.nci.nih.gov/download/attachments/160989212/Div3_2DStructures_Oct2014.sdf']
    filenames = ['structures.sdf']

if __name__ == '__main__':
    DiversityDownloader.download(os.getcwd())
