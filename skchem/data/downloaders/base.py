#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os

from fuel.downloaders.base import default_downloader

class Downloader(object):

    urls = []
    filenames = []

    @classmethod
    def fill_subparser(cls, subparser):
        subparser.set_defaults(urls=cls.urls, filenames=cls.filenames)
        return default_downloader

    @classmethod
    def download(cls, directory=None):
        if not directory:
            directory = os.getcwd()
        return default_downloader(directory, cls.urls, cls.filenames)
