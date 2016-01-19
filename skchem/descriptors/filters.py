#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.descriptors.filters

Chemical filters are defined.

"""

organic = ['H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I']

def is_organic(m):
    return all(atom.element in organic for atom in m.atoms)
