#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.io.sdf

Defining input and output operations for sdf files.
"""

from functools import wraps
import warnings

from rdkit import Chem
import pandas as pd

from ..core import Mol
from ..utils import Suppressor, squeeze

def _drop_props(row):
    for prop in row.structure.props.keys():
        row.structure.ClearProp(prop)

def _set_props(row, cols):
    for i in cols:
        row.structure.SetProp(str(i), str(row[i])) # rdkit props can only be str

def _set_name(row):
    row.structure.name = str(row.name) # rdkit props can only be strs

def read_sdf(sdf, error_bad_mol=False, warn_bad_mol=True, nmols=None,
             skipmols=None, skipfooter=None, read_props=True, mol_props=False,
             *args, **kwargs):

    """Read an sdf file into a `pd.DataFrame`.

    The function wraps the RDKit `ForwardSDMolSupplier` object.

    Args:
        sdf (str or file-like):
            The location of data to load, as a file path, or a file-like object.
        error_bad_mol (bool):
            Whether an error should be raised if a molecule fails to parse.
            Default is False.
        warn_bad_mol (bool):
            Whether a warning should be output if a molecule fails to parse.
            Default is True.
        nmols (int):
            The number of molecules to read. If `None`, read all molecules.
            Default is `None`.
        skipmols (int):
            The number of molecules to skip at start.
            Default is `0`.
        skipfooter (int):
            The number of molecules to skip from the end.
            Default is `0`.
        read_props (bool):
            Whether to read the properties into the data frame.
            Default is `True`.
        mol_props (bool):
            Whether to keep properties in the molecule dictionary after they are
            extracted to the dataframe.
            Default is `False`.
        args, kwargs:
            Arguments will be passed to rdkit's ForwardSDMolSupplier.

    Returns:
        pandas.DataFrame:
            The loaded data frame, with Mols supplied in the `structure` field.

    See also:
        rdkit.Chem.SDForwardMolSupplier
        skchem.read_smiles
    """

    # nmols is actually the index to cutoff.  If we skip some at start, we need
    # to add this number
    if skipmols:
        nmols += skipmols

    if isinstance(sdf, str):
        sdf = open(sdf, 'rb') # use read bytes for python 3 compatibility

    # use the suppression context manager to not pollute our stdout with rdkit
    # errors and warnings.
    # perhaps this should be captured better by Mol etc.
    with Suppressor():

        mol_supp = Chem.ForwardSDMolSupplier(sdf, *args, **kwargs)

        mols = []

        # single loop through sdf
        for i, mol in enumerate(mol_supp):

            if skipmols and i < skipmols:
                continue

            if nmols and i >= nmols:
                break

            # rdkit returns None if it fails to parse a molecule.  We will raise
            # errors unless force is used.
            if mol is None:
                msg = 'Molecule {} could not be decoded.'.format(i + 1)
                if error_bad_mol:
                    raise ValueError(msg)
                elif warn_bad_mol:
                    warnings.warn(msg)
                continue

            mols.append(Mol(mol))


        if skipfooter:
            mols = mols[:-skipfooter]

    idx = pd.Index((m.name for m in mols), name='name')
    data = pd.DataFrame(mols, columns=['structure'])

    if read_props:
        props = pd.DataFrame([{k: v for (k, v) in mol.props.items()} for mol in mols])
        data = pd.concat([data, props], axis=1)
        # now we have extracted the props, we can delete if required
        if not mol_props:
            data.apply(_drop_props, axis=1)

    data.index = idx
    return squeeze(data, axis=1)

def write_sdf(data, sdf, write_cols=True, index_as_name=True, mol_props=False,
              *args, **kwargs):

    """ Write an sdf file from a dataframe.

    Args:
        data (pandas.Series or pandas.DataFrame):
            Pandas data structure with a `structure` column containing compounds
            to serialize.
        sdf (str or file-like):
            A file path or file-like object specifying where to write the
            compound data.
        write_cols (bool):
            Whether columns should be written as props. Default `True`.
        index_as_name (bool):
            Whether to use index as the header, or the molecule's name.
            Default is `True`.
        mol_props (bool):
            Whether to write properties in the Mol dictionary in addition to
            fields in the frame.

    Warn:
        This function will change the names of the compounds if the
        `index_as_name` argument is `True`, and will delete all properties in
        the molecule dictionary if `mol_props` is `False`.
    """
    if isinstance(data, pd.Series):
        data = data.to_frame(name='structure')

    names = [m.name for m in data.structure]

    writer = Chem.SDWriter(sdf, *args, **kwargs)

    cols = list(data.columns.drop('structure'))

    if not mol_props:
        data.apply(_drop_props, axis=1)

    if write_cols:
        data.apply(_set_props, cols=cols, axis=1)

    if index_as_name:
        data.apply(_set_name, axis=1)

    data.structure.apply(writer.write)

    # rdkit writer changes names sometimes
    for mol, name in zip(data.structure, names):
        mol.name = name

@wraps(write_sdf)
def _to_sdf_series(self, *args, **kwargs):

    return write_sdf(self, write_cols=False, *args, **kwargs)

@wraps(write_sdf)
def _to_sdf_df(self, *args, **kwargs):

    return write_sdf(self, *args, **kwargs)

pd.Series.to_sdf = _to_sdf_series
pd.DataFrame.to_sdf = _to_sdf_df


@classmethod
@wraps(read_sdf)
def _from_sdf_df(_, *args, **kwargs):

    return read_sdf(*args, **kwargs)

pd.DataFrame.from_sdf = _from_sdf_df

@classmethod
@wraps(read_sdf)
def _from_sdf_series(_, *args, **kwargs):

    return read_sdf(*args, **kwargs).structure

pd.Series.from_sdf = _from_sdf_series
