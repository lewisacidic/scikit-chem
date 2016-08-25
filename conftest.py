#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# conftest

Configuring pytest for scikit-chem.
"""

def pytest_addoption(parser):
    parser.addoption('--with-chemaxon', action='store_true',
                     help='mark tests that fail if no chemaxon provided.')

    parser.addoption('--run-slow', action='store_true',
                     help='mark tests that take a while to run.')