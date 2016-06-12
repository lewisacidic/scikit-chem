#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class BradleyOpenMPDownloader(Downloader):
    urls = ['https://ndownloader.figshare.com/files/1503990']
    filenames = ['bradley_melting_point_dataset.xlsx']

if __name__ == '__main__':
    BradleyOpenMPDownloader.download(os.getcwd())
