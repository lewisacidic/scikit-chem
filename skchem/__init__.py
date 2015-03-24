#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" A cheminformatics library aiming to integrate into the Scientific Python Stack """

__version__ = '0.0.1dev'

from .core import Mol
from .io import read_sdf, read_smiles
