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

def resource(*args):

    """ passes a file path for a data resource specified """

    return os.path.join(os.path.dirname(__file__), *args)

PERIODIC_TABLE = pd.read_csv(resource('atom_data.csv'), index_col=0)

# must go below periodic table
from .datasets import (
    BursiAmes,
    MullerAmes,
    PhysProp,
    BradleyOpenMP,
    NMRShiftDB2,
    Tox21
)