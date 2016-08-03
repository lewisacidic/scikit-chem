#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
# skchem.base

Base classes for scikit-chem objects.
"""
import subprocess
from abc import ABCMeta, abstractmethod
from tempfile import NamedTemporaryFile
import time
import logging

import pandas as pd

from .utils import NamedProgressBar
from . import core
from .utils import iterable_to_series, optional_second_method, nanarray, squeeze
from . import io

LOGGER = logging.getLogger(__name__)


class BaseTransformer(metaclass=ABCMeta):

    def __init__(self, verbose=True):
        self.verbose = verbose


    def optional_bar(self, **kwargs):
        if self.verbose:
            bar = NamedProgressBar(name=self.__class__.__name__, **kwargs)
        else:
            def bar(x):
                return x
        return bar

    @property
    @abstractmethod
    def axes_names(self):
        pass

    @abstractmethod
    def transform(self, mols):
        pass


class Transformer(BaseTransformer, metaclass=ABCMeta):

    """ Molecular based Transformer object. """

    @property
    @abstractmethod
    def columns(self):
        return pd.Index(None)

    @abstractmethod
    def _transform_mol(self, mol):
        """ Transform a molecule. """
        pass

    def _transform_series(self, ser):
        """ Transform a series of molecules to a list. """

        bar = self.optional_bar()

        return [self._transform_mol(mol) for mol in bar(ser)]

    @optional_second_method
    def transform(self, mols, **kwargs):
        if isinstance(mols, core.Mol):
            # just squeeze works on series
            return pd.Series(self._transform_mol(mols),
                             index=self.columns,
                             name=self.__class__.__name__).squeeze()

        elif not isinstance(mols, pd.Series):
            mols = iterable_to_series(mols)

        res = pd.DataFrame(self._transform_series(mols),
                           index=mols.index,
                           columns=self.columns)

        return squeeze(res, axis=1)

    @property
    def axes_names(self):
        return 'batch', self.columns.name


class BatchTransformer(BaseTransformer, metaclass=ABCMeta):
    """ Transformer in which transforms on multiple molecules save overhead.

    Override `_transform_series` with the transformation rather than `_transform_mol`."""

    def _transform_mol(self, mol):
        v = self.verbose
        self.verbose = False
        res = self.transform([mol]).iloc[0]
        self.verbose = v
        return res

    @abstractmethod
    def _transform_series(self, ser):
        pass


class AtomTransformer(BaseTransformer, metaclass=ABCMeta):
    """ Transformer that will produce a Panel. """

    def __init__(self, max_atoms=100, **kwargs):
        self.max_atoms = max_atoms
        self.major_axis = pd.RangeIndex(self.max_atoms, name='atom_idx')
        super(AtomTransformer, self).__init__(**kwargs)

    @property
    @abstractmethod
    def minor_axis(self):

        return pd.Index(None)  # expects a length

    @property
    def axes_names(self):
        return 'batch', 'atom_idx', self.minor_axis.name

    @optional_second_method
    def transform(self, mols):

        if isinstance(mols, core.Atom):
            # just squeeze works on series
            return pd.Series(self._transform_atom(mols),
                             index=self.minor_axis).squeeze()

        elif isinstance(mols, core.Mol):
            res = pd.DataFrame(self._transform_mol(mols),
                               index=self.major_axis[:len(mols.atoms)],
                               columns=self.minor_axis)
            return squeeze(res, axis=1)

        elif not isinstance(mols, pd.Series):
            mols = iterable_to_series(mols)

        res = pd.Panel(self._transform_series(mols),
                       items=mols.index,
                       major_axis=self.major_axis,
                       minor_axis=self.minor_axis)

        return squeeze(res, axis=(1, 2))

    @abstractmethod
    def _transform_atom(self, atom):
        """ Transform an atom to a 1D array of length `len(self.columns)`. """

        pass

    def _transform_mol(self, mol):
        """ Transform a Mol to a 2D array. """

        res = nanarray((len(mol.atoms), len(self.minor_axis)))
        for i, atom in enumerate(mol.atoms):
            res[i] = self._transform_atom(atom)
        return res

    def _transform_series(self, ser):
        """ Transform a Series<Mol> to a 3D array. """

        if self.verbose:
            bar = NamedProgressBar(name=self.__class__.__name__)
        else:
            # use identity.
            def bar(obj):
                return obj

        res = nanarray((len(ser), self.max_atoms, len(self.minor_axis)))
        for i, mol in enumerate(bar(ser)):
            res[i, :len(mol.atoms), :len(self.minor_axis)] = self._transform_mol(mol)
        return res


class External(metaclass=ABCMeta):
    """ Mixin for wrappers of external CLI tools. """

    install_hint = "" # give an explanation of how to install external tool here.

    def __init__(self, **kwargs):
        assert self.validated, 'External tool not installed. ' + self.install_hint
        super(External, self).__init__(**kwargs)

    @property
    def validated(self):
        if not hasattr(self.__class__, '_validated'):
            self.__class__._validated = self.validate_install()
        return self.__class__._validated

    @staticmethod
    @abstractmethod
    def validate_install():
        """ Determine if the external tool is available. """
        pass


class CLIWrapper(External, BaseTransformer, metaclass=ABCMeta):
    """ CLI wrapper. """

    def __init__(self, error_on_fail=False, warn_on_fail=True, **kwargs):
        super(CLIWrapper, self).__init__(**kwargs)
        self.error_on_fail = error_on_fail
        self.warn_on_fail = warn_on_fail

    def _transform_series(self, ser):
        with NamedTemporaryFile(suffix='.sdf') as infile, NamedTemporaryFile() as outfile:
            io.write_sdf(ser, infile.name)
            args = self._cli_args(infile.name, outfile.name)
            p = subprocess.Popen(args, stderr=subprocess.PIPE)

            if self.verbose:
                bar = self.optional_bar(max_value=len(ser))
                while p.poll() is None:
                    time.sleep(0.5)
                    bar.update(self.monitor_progress(outfile.name))
                bar.finish()

            p.wait()
            res = self._parse_outfile(outfile.name)

        errs = p.stderr.read().decode()
        errs = self._parse_errors(errs)
        # set the index of results to that of the input, with the failed indices removed
        if isinstance(res, (pd.Series, pd.DataFrame)):
            res.index = ser.index.delete(errs)
        elif isinstance(res, pd.Panel):
            res.items = ser.index.delete(errs)
        else:
            raise ValueError('Parsed datatype ({}) not supported.'.format(type(res)))

        # go through the errors and put them back in (transform doesn't lose instances)
        if len(errs):
            for err in errs:
                err = ser.index[err]
                if self.error_on_fail:
                    raise ValueError('Failed to transform {}.'.format(err))
                if self.warn_on_fail:
                    LOGGER.warn('Failed to transform %s', err)
                res.ix[err] = None

        return res.loc[ser.index].values

    @abstractmethod
    def _cli_args(self, infile, outfile):
        return []

    @abstractmethod
    def monitor_progress(self, filename):
        """ Report the progress. """
        pass

    @abstractmethod
    def _parse_outfile(self, outfile):
        pass

    @abstractmethod
    def _parse_errors(self, errs):
        pass


class Pipeline(object):

    """ Pipeline object. Applies filters and transformers. """

    def transform(self):
        raise NotImplemented


class Featurizer(metaclass=ABCMeta):

    """ Base class for m -> data transforms, such as Fingerprinting etc."""

    pass

