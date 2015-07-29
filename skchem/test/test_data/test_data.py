#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" Tests for data functions """

from skchem.data import resource

def test_resource():

    """ Does resource target the the methane smiles test file?"""

    with open(resource('methane.smiles'), 'r') as f:
        assert f.readline() == 'C'
