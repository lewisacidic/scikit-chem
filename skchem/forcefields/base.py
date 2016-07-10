#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.base

Module specifying base class for forcefields.
"""

class ForceField(object):
    def __init__(self, warn_on_fail=True, error_on_fail=True):
        self.warn_on_fail = warn_on_fail
        self.error_on_fail = error_on_fail


    def optimize(self, mol):
        raise NotImplementedError