#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields

Module specifying forcefields.
"""

from .mmff import MMFF
from .uff import UFF
from .base import RoughEmbedding

def get(name):
    DEFAULTS = {
        'uff': UFF,
        'mmff': MMFF,
        'rough': RoughEmbedding,
    }
    if isinstance(name, str):
        return DEFAULTS[name]()