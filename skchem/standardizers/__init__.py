#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from .chemaxon import ChemAxonStandardizer

__all__ = [
    'ChemAxonStandardizer'
]


def get(identifier):
    if isinstance(identifier, str):
        defaults = {'chemaxon': ChemAxonStandardizer}
        return defaults[identifier]()
    else:
        return identifier
