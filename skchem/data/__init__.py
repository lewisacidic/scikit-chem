#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.data

Module for handling data. Data can be accessed using the resource function.

"""

import os
import pandas as pd

# must go below periodic table
from .datasets import (
    Diversity,
    BursiAmes,
    MullerAmes,
    PhysProp,
    BradleyOpenMP,
    NMRShiftDB2,
    Tox21
)