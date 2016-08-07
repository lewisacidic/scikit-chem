#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# file title

Description
"""


import logging

from .base import Dataset
from ..downloaders.diversity import DiversityDownloader
from ..converters.diversity_set import DiversityConverter


class Diversity(Dataset):

    """ Example dataset, the NCI DTP Diversity Set III. """

    filename = 'diversity.h5'
    downloader = DiversityDownloader
    converter = DiversityConverter

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    Example.download()

