#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.descriptors.physicochemical

Physicochemical descriptors and associated functions are defined.

"""

from rdkit.Chem.Descriptors import MolWt

def molecular_weight(m):
    return MolWt(m)
