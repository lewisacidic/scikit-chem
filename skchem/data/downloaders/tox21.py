#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import os
from .base import Downloader

class Tox21Downloader(Downloader):
    filenames = ['data_allsdf', 'challenge_testsdf',
                 'challenge_scoresdf', 'challenge_scoretxt']
    urls = ['https://tripod.nih.gov/tox21/challenge/download?id=tox21_10k_' \
            + name for name in filenames]
    filenames = ['train.sdf.zip', 'valid.sdf.zip', 'test.sdf.zip', 'test.txt']

if __name__ == '__main__':
    Tox21Downloader.download(os.getcwd())
