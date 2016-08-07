#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from .bursi_ames import BursiAmesDownloader
from .muller_ames import MullerAmesDownloader
from .tox21 import Tox21Downloader
from .nmrshiftdb2 import NMRShiftDB2Downloader
from .physprop import PhysPropDownloader
from .bradley_open_mp import BradleyOpenMPDownloader
from .diversity import DiversityDownloader

all_downloaders = (
    ('diversity', DiversityDownloader.fill_subparser),
    ('bursi_ames', BursiAmesDownloader.fill_subparser),
    ('muller_ames', MullerAmesDownloader.fill_subparser),
    ('tox21', Tox21Downloader.fill_subparser),
    ('nmrshiftdb2', NMRShiftDB2Downloader.fill_subparser),
    ('physprop', PhysPropDownloader.fill_subparser),
    ('bradley_open_mp', BradleyOpenMPDownloader.fill_subparser)
)
