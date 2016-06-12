#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from .bursi_ames import BursiAmesDownloader
from .muller_ames import MullerAmesDownloader
from .tox21 import Tox21Downloader
#from .nmrshiftdb2 import NMRShiftDB2Downloader
from .physprop import PhysPropDownloader
#from .wombat import WombatDownloader
#from .chembl import ChEMBLDownloader
from .bradley_open_mp import BradleyOpenMPDownloader

__version__ = '0.0.5'

all_downloaders = (
    ('bursi_ames', BursiAmesDownloader.fill_subparser),
    ('muller_ames', MullerAmesDownloader.fill_subparser),
    ('tox21', Tox21Downloader.fill_subparser),
#    ('nmrshiftdb2', NMRShiftDB2Downloader.fill_subparser),
    ('physprop', PhysPropDownloader.fill_subparser),
    ('bradley_open_mp', BradleyOpenMPDownloader.fill_subparser),
#    ('wombat', WombatDownloader.fill_subparser), ## not yet available
#    ('chembl', ChEMBLDownloader.fill_subparser) ## not yet available
)
