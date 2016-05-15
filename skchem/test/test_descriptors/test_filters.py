#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" Tests for skchem.descriptors.filters """

from skchem.data import resource
import skchem as skc
from skchem import filters

class TestIsOrganic(object):

    HYDROCARBON = "CCC"
    CISPLATIN = "N[Pt](N)(Cl)Cl"
    ACICLOVIR = "NC1=NC(=O)C2=C(N1)N(COCCO)C=N2"

    def test_hydrocarbon(self):
        assert filters.is_organic(skc.Mol.from_smiles(self.HYDROCARBON))

    def test_big_drug(self):
        assert filters.is_organic(skc.Mol.from_smiles(self.ACICLOVIR))

    def test_inorganic_drug(self):
        assert not filters.is_organic(skc.Mol.from_smiles(self.CISPLATIN))
