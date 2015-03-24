#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.data

Module for handling data. Data can be accessed using the resource function.

"""

import os

def resource(*args):

    """ passes a file path for a data resource specified """

    return os.path.join(os.path.dirname(__file__), *args)
