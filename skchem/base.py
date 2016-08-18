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
import multiprocessing
from tempfile import NamedTemporaryFile
import time
import logging

import pandas as pd

from .utils import NamedProgressBar, DummyProgressBar
from . import core
from .utils import (iterable_to_series, optional_second_method, nanarray,
                    squeeze, yaml_dump, json_dump)
from . import io

LOGGER = logging.getLogger(__name__)


class BaseTransformer(object):

    """ Transformer Base Class.

    Specific Base Transformer classes inherit from this class and implement
    `transform` and `axis_names`.
    """

    __metaclass__ = ABCMeta

    # To share some functionality betweeen Transformer and AtomTransformer

    def __init__(self, n_jobs=1, verbose=True):
        self._n_jobs = None  # property cache
        self.n_jobs = n_jobs
        self.verbose = verbose

    @property
    def n_jobs(self):
        return self._n_jobs

    @n_jobs.setter
    def n_jobs(self, val):
        if val >= 1:
            self._n_jobs = val
        elif val == -1:
            self._n_jobs = multiprocessing.cpu_count()

    def get_params(self):
        """ Get a dictionary of the parameters of this object. """
        params = list(self.__class__.__init__.__code__.co_varnames)
        params.remove('self')
        return {param: getattr(self, param) for param in params}

    @classmethod
    def from_params(cls, params):
        """ Create a instance from a params dictionary. """
        return cls(**params)

    def to_dict(self):

        """ Return a dictionary representation of the object."""

        n = '{}.{}'.format(self.__class__.__module__, self.__class__.__name__)
        return {n: self.get_params()}

    def to_json(self, target=None):

        """ Serialize the object as JSON.

        Args:
            target (str or file-like):
                A file or filepath to serialize the object to.  If `None`,
                return the JSON as a string.

            Returns:
                None or str
        """

        return json_dump(self.to_dict(), target)

    def to_yaml(self, target=None):

        """ Serialize the object as YAML.

        Args:
            target (str or file-like):
                A file or filepath to serialize the object to.  If `None`,
                return the YAML as a string.

            Returns:
                None or str
        """

        return yaml_dump(self.to_dict(), target)

    def copy(self):
        """ Return a copy of this object. """
        return self.__class__(**self.get_params())

    def optional_bar(self, **kwargs):
        if self.verbose:
            bar = NamedProgressBar(name=self.__class__.__name__, **kwargs)
        else:
            bar = DummyProgressBar(**kwargs)
        return bar

    @property
    @abstractmethod
    def axes_names(self):
        """ tuple: The names of the axes. """
        pass

    @abstractmethod
    def transform(self, mols):
        """ Transform objects according to the objects transform protocol.

        Args:
            mols (skchem.Mol or pd.Series or iterable):
                The mol objects to transform.

        Returns:
            pd.Series or pd.DataFrame
        """
        pass

    def __eq__(self, other):
        return self.get_params() == other.get_params()


class Transformer(BaseTransformer):

    """ Molecular based Transformer Base class.

    Concrete Transformers inherit from this class and must implement
    `_transform_mol` and `_columns`.

    See Also:
         AtomTransformer."""

    @property
    @abstractmethod
    def columns(self):
        """ pd.Index: The column index to use. """
        return pd.Index(None)

    @abstractmethod
    def _transform_mol(self, mol):
        """ Transform a molecule. """
        pass

    def _transform_series(self, ser):
        """ Transform a series of molecules to an np.ndarray. """
        LOGGER.debug('Transforming series of length %s with %s jobs',
                     len(ser), self.n_jobs)

        bar = self.optional_bar(max_value=len(ser))
        if self.n_jobs == 1:
            return [self._transform_mol(mol) for mol in bar(ser)]
        else:
            cpy = self.copy()
            with multiprocessing.Pool(processes=self.n_jobs) as pool:
                return [res for res in bar(pool.imap(cpy._transform_mol, ser))]

    @optional_second_method
    def transform(self, mols, **kwargs):
        """ Transform objects according to the objects transform protocol.

        Args:
            mols (skchem.Mol or pd.Series or iterable):
                The mol objects to transform.

        Returns:
            pd.Series or pd.DataFrame
        """
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
        """ tuple: The names of the axes. """
        return 'batch', self.columns.name


class BatchTransformer(BaseTransformer):
    """  Mixin for which transforms on multiple molecules save overhead.

    Implement `_transform_series` with the transformation rather than
    `_transform_mol`. Must occur before `Transformer` or  `AtomTransformer` in
    method resolution order.

    See Also:
         Transformer, AtomTransformer.
    """

    def _transform_mol(self, mol):
        """ Transform a molecule. """

        v = self.verbose
        self.verbose = False
        res = self.transform([mol]).iloc[0]
        self.verbose = v
        return res

    @abstractmethod
    def _transform_series(self, ser):
        """ Transform a series of molecules to an np.ndarray. """
        pass


class AtomTransformer(BaseTransformer):
    """ Transformer that will produce a Panel.

    Concrete classes inheriting from this should implement `_transform_atom`,
    `_transform_mol` and `minor_axis`.

    See Also:
        Transformer
    """

    def __init__(self, max_atoms=100, **kwargs):
        self.max_atoms = max_atoms
        self.major_axis = pd.RangeIndex(self.max_atoms, name='atom_idx')
        super(AtomTransformer, self).__init__(**kwargs)

    @property
    @abstractmethod
    def minor_axis(self):
        """ pd.Index: Minor axis of transformed values.  """
        return pd.Index(None)  # expects a length

    @property
    def axes_names(self):
        """ tuple: The names of the axes. """
        return 'batch', 'atom_idx', self.minor_axis.name

    @optional_second_method
    def transform(self, mols):
        """ Transform objects according to the objects transform protocol.

        Args:
            mols (skchem.Mol or pd.Series or iterable):
                The mol objects to transform.

        Returns:
            pd.Series or pd.DataFrame
        """
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
        LOGGER.debug('Transforming series of length %s with %s jobs',
                     len(ser), self.n_jobs)
        bar = self.optional_bar(max_value=len(ser))

        res = nanarray((len(ser), self.max_atoms, len(self.minor_axis)))

        if self.n_jobs == 1:
            for i, mol in enumerate(bar(ser)):
                res[i, :len(mol.atoms),
                    :len(self.minor_axis)] = self._transform_mol(mol)
        else:
            cpy = self.copy()
            with multiprocessing.Pool(self.n_jobs) as pool:
                for (i, ans) in enumerate(bar(pool.imap(cpy._transform_mol,
                                                        ser))):
                    res[i, :len(ans), :len(self.minor_axis)] = ans
        return res

class External(object):
    """ Mixin for wrappers of external CLI tools.

    Concrete classes must implement `validate_install`.

    Attributes:
        install_hint (str): an explanation of how to install external tool.
    """

    __metaclass__ = ABCMeta

    install_hint = ""

    def __init__(self, **kwargs):
        if not self.validated:
            msg = 'External tool not installed. {}'.format(self.install_hint)
            raise RuntimeError(msg)
        super(External, self).__init__(**kwargs)

    @property
    def validated(self):
        """ bool: whether the external tool is installed and active. """
        if not hasattr(self.__class__, '_validated'):
            self.__class__._validated = self.validate_install()
        return self.__class__._validated

    @staticmethod
    @abstractmethod
    def validate_install():
        """ Determine if the external tool is available. """
        pass


class CLIWrapper(External, BaseTransformer):
    """ CLI wrapper.

    Concrete classes inheriting from this must implement `_cli_args`,
    `monitor_progress`, `_parse_outfile`, `_parse_errors`."""

    def __init__(self, error_on_fail=False, warn_on_fail=True, **kwargs):
        super(CLIWrapper, self).__init__(**kwargs)
        self.error_on_fail = error_on_fail
        self.warn_on_fail = warn_on_fail

    @property
    def n_jobs(self):
        return self._n_jobs

    @n_jobs.setter
    def n_jobs(self, val):
        if val != 1:
            raise NotImplementedError('Multiprocessed external code is not yet'
                                      ' supported.')
        else:
            self._n_jobs = val

    def _transform_series(self, ser):
        """ Transform a series. """
        with NamedTemporaryFile(suffix='.sdf') as infile, \
                NamedTemporaryFile() as outfile:
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
        # set the index of results to that of the input, with the failed
        # indices removed
        if isinstance(res, (pd.Series, pd.DataFrame)):
            res.index = ser.index.delete(errs)
        elif isinstance(res, pd.Panel):
            res.items = ser.index.delete(errs)
        else:
            msg = 'Parsed datatype ({}) not supported.'.format(type(res))
            raise ValueError(msg)

        # go through the errors and put them back in
        # (transform doesn't lose instances)
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
        """ list: The cli arguments. """
        return []

    @abstractmethod
    def monitor_progress(self, filename):
        """ Report the progress. """
        pass

    @abstractmethod
    def _parse_outfile(self, outfile):
        """ Parse the file written and return a series. """
        pass

    @abstractmethod
    def _parse_errors(self, errs):
        """ Parse stderr and return error indices. """
        pass


class Featurizer(object):

    """ Base class for m -> data transforms, such as Fingerprinting etc.

    Concrete subclasses should implement `name`, returning a string uniquely
    identifying the featurizer. """

    __metaclass__ = ABCMeta
