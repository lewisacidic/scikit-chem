#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.standardizers.chemaxon

Module wrapping ChemAxon Standardizer.  Must have standardizer installed and
license activated.
"""

import os
from tempfile import NamedTemporaryFile
import subprocess
import logging

logger = logging.getLogger(__file__)

import pandas as pd

from .. import core
from .. import io

# ideally we will programatically build this file, but for now just use it.
DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), 'default_config.xml')

class ChemAxonStandardizer(object):

    """ Object wrapping the ChemAxon Standardizer, for standardizing molecules.

    Args:
        config_path (str):
            The path of the config_file. If None, use the default one.
    """
    def __init__(self, config_path=None, keep_failed=True):

        if not config_path:
            config_path = DEFAULT_CONFIG
        self.config_path = config_path
        self.keep_failed = keep_failed

    def transform(self, obj):
        if isinstance(obj, core.Mol):
            return self._transform_mol(obj)
        elif isinstance(obj, pd.Series):
            return self._transform_ser(obj)
        elif isinstance(obj, pd.DataFrame):
            obj = obj.copy()
            res = self._transform_ser(obj.structure)
            obj.iloc['structure'] = res
            return obj

        else:
            raise NotImplementedError

    def _transform_mol(self, mol):
        mol = pd.DataFrame([mol], index=[mol.name], columns=['structure'])
        return self.transform(mol).structure.iloc[0]


    def _transform_by(self, X, by='sdf'):

        with NamedTemporaryFile() as f_in, NamedTemporaryFile() as f_out:
            getattr(io, 'write_' + by)(X, f_in.name)
            args = ['standardize', f_in.name,
                             '-c', self.config_path,
                             '-f', 'sdf',
                             '-o', f_out.name,
                             '--ignore-error']
            logger.debug('Running %s', ' '.join(args))
            subprocess.call(args)
            out = getattr(io, 'read_' + by)(f_out.name).structure
        return out


    def _transform_ser(self, X, y=None):
        # TODO: should check for errors, and schedule multiple runs with
        # different serialization methods here.

        out = self._transform_by(X)
        return out

        # failed = X.ix[failed]
        # now_succ = self._transform_by(failed, 'smiles')
        # still_failed = set(now_succ.index).difference(set(failed.index))
        #
        # for fail in still_failed:
        #     logger.warn('%s failed to standardize.', fail)
        #
        # out = out.append(now_succ)
        # if self.keep_failed:
        #     out = out.append(X.ix[still_failed])
        # return out.structure
