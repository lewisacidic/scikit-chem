#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from skchem.tests.utils import resource

def test_resource():

    """ Does resource target the the methane smiles test file?"""

    with open(resource('smiles', 'single_molecule.smiles'), 'r') as f: 
        assert str(f.readline()) == 'C'
