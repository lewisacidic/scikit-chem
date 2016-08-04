#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.io.smiles

Defining input and output operations for smiles files.
"""

import warnings
from functools import wraps

import pandas as pd

from ..utils import Suppressor
from ..core import Mol

def read_smiles(smiles_file, smiles_column=0, name_column=None, delimiter='\t',
                title_line=False, error_bad_mol=False, warn_bad_mol=True,
                drop_bad_mol=True, *args, **kwargs):

    """Read a smiles file into a pandas dataframe.

    The class wraps the pandas read_csv function.

    smiles_file (str, file-like):
        Location of data to load, specified as a string or passed directly as a
        file-like object.  URLs may also be used, see the pandas.read_csv
        documentation.
    smiles_column (int):
        The column index at which SMILES are provided.
        Defaults to `0`.
    name_column (int):
        The column index at which compound names are provided, for use as the
        index in the DataFrame.  If None, use the default index.
        Defaults to `None`.
    delimiter (str):
        The delimiter used.
        Defaults to `\\t`.
    title_line (bool):
        Whether a title line is provided, to use as column titles.
        Defaults to `False`.
    error_bad_mol (bool):
        Whether an error should be raised when a molecule fails to parse.
        Defaults to `False`.
    warn_bad_mol (bool):
        Whether a warning should be raised when a molecule fails to parse.
        Defaults to `True`.
    drop_bad_mol (bool):
        If true, drop any column with smiles that failed to parse. Otherwise,
        the field is None. Defaults to `True`.
    args, kwargs:
        Arguments will be passed to pandas read_csv arguments.

    Returns:
        pandas.DataFrame:
            The loaded data frame, with Mols supplied in the `structure` field.

    See Also:
        pandas.read_csv
        skchem.Mol.from_smiles
        skchem.io.sdf
    """

    with Suppressor():

        # set the header line to pass to the pandas parser
        # we accept True as being line zero, as is usual for smiles
        # if user specifies a header already, then do nothing

        header = kwargs.pop('header', None)
        if title_line is True:
            header = 0
        elif header is not None:
            pass #remove from the kwargs to not pass it twice
        else:
            header = None

        # read the smiles file
        data = pd.read_csv(smiles_file, delimiter=delimiter, header=header,
                           *args, **kwargs)

        # replace the smiles column with the structure column
        lst = list(data.columns)
        lst[smiles_column] = 'structure'
        data.columns = lst

        def parse(row):
            """ Parse smiles for row """
            try:
                return Mol.from_smiles(row.structure)
            except ValueError:
                msg = 'Molecule {} could not be decoded.'.format(row.name)
                if error_bad_mol:
                    raise ValueError(msg)
                elif warn_bad_mol:
                    warnings.warn(msg)

                return None

        data['structure'] = data['structure'].apply(str)
        data['structure'] = data.apply(parse, axis=1)

        if drop_bad_mol:
            data = data[data['structure'].notnull()]

        # set index if passed
        if name_column is not None:
            data = data.set_index(data.columns[name_column])

        return data


def write_smiles(data, smiles_path):

    """ Write a dataframe to a smiles file.

    Args:
        data (pd.Series or pd.DataFrame):
            The dataframe to write.
        smiles_path (str):
            The path to write the dataframe to.
    """

    if isinstance(data, pd.Series):
        data = data.to_frame(name='structure')
    data['structure'] = data.structure.apply(lambda m: m.to_smiles())
    data = data.reset_index()
    cols = list(data.columns)
    cols.insert(0, cols.pop(cols.index('structure')))
    data = data.reindex(columns=cols)[cols]
    data.to_csv(smiles_path, sep='\t', header=None, index=None)


@classmethod
@wraps(read_smiles)
def _from_smiles_df(_, *args, **kwargs):
    return read_smiles(*args, **kwargs)

@classmethod
@wraps(read_smiles)
def _from_smiles_series(_, *args, **kwargs):
    return read_smiles(*args, **kwargs).structure

@wraps(write_smiles)
def _to_smiles_df(self, *args, **kwargs):
    return write_smiles(self, *args, **kwargs)

pd.DataFrame.from_smiles = _from_smiles_df
pd.Series.from_smiles = _from_smiles_series
pd.Series.to_smiles = _to_smiles_df
pd.DataFrame.to_smiles = _to_smiles_df
