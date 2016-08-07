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
import sys
import re
import subprocess
import logging
import warnings

import pandas as pd

from .. import io
from ..utils import sdf_count
from ..base import CLIWrapper, Transformer, BatchTransformer
from ..filters.base import TransformFilter

LOGGER = logging.getLogger(__name__)

if sys.version_info[0] == 2:
    NoFoundError = OSError
    subprocess.DEVNULL = open(os.devnull, 'w')
else:
    NoFoundError = FileNotFoundError


class ChemAxonStandardizer(CLIWrapper, BatchTransformer, Transformer, TransformFilter):

    """ ChemAxon Standardizer Wrapper.

    Args:
        config_path (str):
            The path of the config_file. If None, use the default one.

    Notes:
        ChemAxon Standardizer must be installed and accessible as `standardize`
        from the shell launching the program.

    Warnings:
        Must use a unique index (see #31).

    Examples:

        >>> import skchem
        >>> std = skchem.standardizers.ChemAxonStandardizer() # doctest:+SKIP
        >>> m = skchem.Mol.from_smiles('CC.CCC')
        >>> print(std.transform(m)) # doctest:+SKIP
        <Mol: CCC>

        >>> data = [m, skchem.Mol.from_smiles('C=CO'), skchem.Mol.from_smiles('C[O-]')]
        >>> std.transform(data) # doctest:+SKIP
        0     <Mol: CCC>
        1    <Mol: CC=O>
        2      <Mol: CO>
        Name: structure, dtype: object

        >>> will_fail = mol = '''932-97-8
        ...      RDKit          3D
        ...
        ...   9  9  0  0  0  0  0  0  0  0999 V2000
        ...    -0.9646    0.0000    0.0032 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...    -0.2894   -1.2163    0.0020 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...    -0.2894    1.2163    0.0025 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...    -2.2146    0.0000   -0.0004 N   0  0  0  0  0  0  0  0  0  0  0  0
        ...     1.0710   -1.2610    0.0002 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...     1.0710    1.2610    0.0007 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...    -3.3386    0.0000   -0.0037 N   0  0  0  0  0  0  0  0  0  0  0  0
        ...     1.8248    0.0000   -0.0005 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...     3.0435    0.0000   -0.0026 O   0  0  0  0  0  0  0  0  0  0  0  0
        ...   1  2  1  0
        ...   1  3  1  0
        ...   1  4  2  3
        ...   2  5  2  0
        ...   3  6  2  0
        ...   4  7  2  0
        ...   5  8  1  0
        ...   8  9  2  0
        ...   6  8  1  0
        ... M  CHG  2   4   1   7  -1
        ... M  END
        ... '''

        >>> will_fail = skchem.Mol.from_molblock(will_fail)
        >>> std.transform(will_fail) # doctest:+SKIP
        nan

        >>> data = [will_fail] + data

        >>> std.transform(data) # doctest:+SKIP
        0           None
        1     <Mol: CCC>
        2    <Mol: CC=O>
        3      <Mol: CO>
        Name: structure, dtype: object

        >>> std.transform_filter(data) # doctest:+SKIP
        1     <Mol: CCC>
        2    <Mol: CC=O>
        3      <Mol: CO>
        Name: structure, dtype: object

        >>> std.keep_failed = True # doctest:+SKIP
        >>> std.transform(data) # doctest:+SKIP
        0    <Mol: [N-]=[N+]=C1C=CC(=O)C=C1>
        1                         <Mol: CCC>
        2                        <Mol: CC=O>
        3                          <Mol: CO>
        Name: structure, dtype: object

    """
    install_hint = """ Install ChemAxon from https://www.chemaxon.com.  It requires a license,
    which can be freely obtained for academics. """

    DEFAULT_CONFIG = os.path.join(os.path.dirname(__file__), 'default_config.xml')

    def __init__(self, config_path=None, keep_failed=False, **kwargs):

        super(ChemAxonStandardizer, self).__init__(**kwargs)

        if not config_path:
            config_path = self.DEFAULT_CONFIG
        self.config_path = config_path
        self.keep_failed = keep_failed

    @property
    def columns(self):
        return ['structure']

    def _transform_series(self, ser):

        # implement keep_failed functionality here
        res = super(ChemAxonStandardizer, self)._transform_series(ser)
        mask = pd.isnull(res)

        for m_in, m_out in zip(ser[~mask], res[~mask]):
            m_out.name = m_in.name

        if self.keep_failed:
            res[mask] = ser.iloc[mask]
        return res

    def _parse_outfile(self, outfile):
        """ Reads output file and returns a list"""
        return io.read_sdf(outfile, read_props=False)

    def _parse_errors(self, errs):
        """ Reads stderr and parses out failures as a list of indices. """
        LOGGER.debug('stderr: %s', errs if errs else None)
        errs = errs.strip().split('\n')
        errs = [re.findall('No. ([0-9]+):', err) for err in errs]
        return [int(err[0]) - 1 for err in errs if len(err)]

    def _cli_args(self, infile, outfile):
         """ The command line arguments to use for the subprocess. """

         return  ['standardize', infile,
                '-c', self.config_path,
                '-f', 'sdf',
                '-o', outfile,
                '--ignore-error']

    @staticmethod
    def validate_install():
        """ Check if we can call cxcalc. """
        try:
            return subprocess.call(['standardize', '-h'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
        except NoFoundError:
            return False

    def monitor_progress(self, filename):
        return sdf_count(filename)

    def filter(self, *args, **kwargs):
        warnings.warn('Filter returns the unstandardized Mols. Did you mean to use `transform_filter`?')
        super(ChemAxonStandardizer, self).filter(*args, **kwargs)