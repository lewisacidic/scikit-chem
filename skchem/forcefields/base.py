#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.base

Module specifying base class for forcefields.
"""
import warnings
import pandas as pd
import progressbar
from rdkit.Chem.rdDistGeom import EmbedMolecule

from .. import core

class ForceField(object):
    def __init__(self, embed=True, warn_on_fail=True, error_on_fail=False, drop_failed=True, add_hs=True):

        self.add_hs = add_hs
        self.drop_failed = drop_failed
        self.warn_on_fail = warn_on_fail
        self.error_on_fail = error_on_fail
        self.preembed = embed

    def embed(self, mol):

        success = EmbedMolecule(mol)
        if success == -1:
            msg = 'Failed to Embed Molecule {}'.format(mol.name)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
                return None
            else:
                pass

        if self.add_hs:
            return mol.add_hs(add_coords=True)
        else:
            return mol


    def optimize(self, mol):

        # TODO: likely need to handle which conformer here

        if self.preembed:
            mol = self.embed(mol)

        if mol is None: # embedding failed
            return None

        res = self._optimize(mol)

        if res == -1:
            msg = 'Failed to optimize molecule \'{}\' using {}'.format(mol.name, self.__class__)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
            return None

        return mol

    def _optimize(self, mol):
        raise NotImplementedError

    def transform(self, obj):
        if isinstance(obj, core.Mol):
            return self.optimize(obj)
        elif isinstance(obj, pd.Series):
            bar = progressbar.ProgressBar()
            res = pd.Series([self.optimize(mol) for mol in bar(obj)], index=obj.index, name=obj.name)
            if self.drop_failed:
                res = res[res.notnull()]
            return res
        elif isinstance(obj, pd.DataFrame):
            res = self.transform(obj.structure)
            return obj.ix[res.index]
        elif isinstance(obj, (tuple, list)):
            return self.transform(pd.Series(obj, [mol.name for mol in obj]))
        else:
            raise NotImplementedError

class RoughEmbedding(ForceField):
    def _optimize(self, mol):
        return mol
