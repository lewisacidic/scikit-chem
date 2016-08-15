#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.data.downloaders.chembl

ChEMBL dataset downloader
"""

from .base import Downloader


class ChEMBLDownloader(Downloader):
    urls = []
    filenames = ['chembl_raw.h5']

    def __init__(self):
        raise NotImplementedError
        super(ChEMBLDownloader, self).__init__()

