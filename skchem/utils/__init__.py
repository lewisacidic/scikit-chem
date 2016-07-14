#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.utils

Module providing utility functions for scikit-chem
"""

from .suppress import Suppressor
from .string import camel_to_snail
from .decorators import (
    takes_mol_series,
    method_takes_mol_series,
    takes_pandas,
    method_takes_pandas
)
from .progressbar import NamedProgressBar
from .io import line_count