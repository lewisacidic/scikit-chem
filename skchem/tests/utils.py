#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
Utilities for testing.
"""
import os

def resource(*names):
    return os.path.join(os.path.split(os.path.abspath(__file__))[0], 'resources', *names)
