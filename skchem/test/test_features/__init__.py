#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.test.test_features

Tests for the `features` package.
"""

import pytest
from ... import Mol


def phenethylamine(a_x='H', a_y='H'):

    a_x = '' if a_x == 'H' else '({})'.format(a_x)
    a_y = '' if a_y == 'H' else '({})'.format(a_y)

    return Mol.from_smiles('CN(C)CC(Br)c1cc{}c{}cc1'.format(a_x, a_y))


ATS = [('H', 'H'),
       ('H', 'F'), ('H', 'Cl'), ('H', 'Br'), ('H', 'I'), ('H', 'C'),
       ('F', 'H'), ('Cl', 'H'), ('Br', 'H'), ('I', 'H'), ('C', 'H'),
       ('Cl', 'F'), ('Br', 'F'), ('C', 'F'),
       ('Cl', 'Cl'), ('Br', 'Cl'), ('C', 'Cl'),
       ('Cl', 'Br'), ('Br', 'Br'), ('C', 'Br'),
       ('C', 'C'), ('Br', 'C')]


PHENS = [phenethylamine(x, y) for x, y in ATS]

