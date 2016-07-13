#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.io

Module defining input and output methods in scikit-chem.

"""

from .sdf import read_sdf, write_sdf
from .smiles import read_smiles
