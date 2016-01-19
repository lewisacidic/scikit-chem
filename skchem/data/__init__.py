#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.data

Module for handling data. Data can be accessed using the resource function.

"""

import os
import pandas as pd

def resource(*args):

    """ passes a file path for a data resource specified """

    return os.path.join(os.path.dirname(__file__), *args)

def periodic_table():

    return pd.read_csv(resource('atomic_data.csv'))
