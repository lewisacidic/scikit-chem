#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from .chemaxon import ChemAxonStandardizer


def get(identifier):
    if isinstance(identifier, str):
        DEFAULTS = {'chemaxon': ChemAxonStandardizer}
        return DEFAULTS[identifier]()
    else:
        return identifier
