#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from .bursi_ames import BursiAmesConverter
from .muller_ames import MullerAmesConverter
from .tox21 import Tox21Converter
from .nmrshiftdb2 import NMRShiftDB2Converter
from .physprop import PhysPropConverter
from .bradley_open_mp import BradleyOpenMPConverter
from .diversity_set import DiversityConverter

all_converters = (
    ('diversity', DiversityConverter.fill_subparser),
    ('bursi_ames', BursiAmesConverter.fill_subparser),
    ('muller_ames', MullerAmesConverter.fill_subparser),
    ('tox21', Tox21Converter.fill_subparser),
    ('nmrshiftdb2', NMRShiftDB2Converter.fill_subparser),
    ('physprop', PhysPropConverter.fill_subparser),
    ('bradley_open_mp', BradleyOpenMPConverter.fill_subparser)
)
