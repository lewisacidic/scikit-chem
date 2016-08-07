#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# file title

Description
"""

import pytest

def pytest_addoption(parser):
    parser.addoption("--chemaxon", action="store_true",
        help="mark tests that fail if no chemaxon provided.")
