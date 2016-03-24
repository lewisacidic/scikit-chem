#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.io.sdf

Defining input and output operations for sdf files.
"""

from functools import wraps
import warnings

from rdkit import Chem
import skchem
from skchem.utils import Suppressor
import pandas as pd

def _drop_props(row):
    for prop in row.structure.props.keys():
        row.structure.ClearProp(prop)

def _set_props(row, cols):
    for i in cols:
        row.structure.SetProp(str(i), str(row[i])) # rdkit props can only be strs

def _set_name(row):
    row.structure.name = str(row.name) # rdkit props can only be strs

def read_sdf(sdf, error_bad_mol=False, warn_bad_mol=True, nmols=None,
             skipmols=None, skipfooter=None, read_props=True, mol_props=False,
             *args, **kwargs):

    """
        Read an sdf file into a pandas dataframe.
        The function wraps the RDKit ForwardSDMolSupplier object.

        @param sdf           A file path provided as a :str:, or a :file-like:
                             object.
        @param error_bad_mol A :bool: specifying if an error should be raised if
                             a molecule fails to parse.
        @param warn_bad_mol  A :bool: specifying if a warning should be output
                             if a molecule fails to parse.
        @param nmols         An :int: specifying number of molecules to read.
                             If none, read all molecules.
        @param skipmols      An :int: specifying number of molecules to skip at
                             start.
        @param skipfooter    An :int: specifying number of molecules to skip
                             from the end.
        @param mol_props     A :bool: specifying whether to keep properties in
                             the molecule dictionary.
        Additionally, ForwardSDMolSupplier arguments will be passed.

        @returns df         A dataframe of type :pandas.core.frame.DataFrame:.

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

            mols.append(skchem.Mol(mol))


        if skipfooter:
            mols = mols[:-skipfooter]

    idx = pd.Index((m.name for m in mols), name='name')
    data = pd.DataFrame(mols, columns=['structure'])

    if read_props:
        props = pd.DataFrame([mol.props for mol in mols])
        data = pd.concat([data, props], axis=1)

    # now we have extracted the props, we can delete if required
    if not mol_props:
        data.apply(_drop_props, axis=1)

    data.index = idx
    return data

def write_sdf(df, sdf, write_cols=True, index_as_name=True, mol_props=False,
              *args, **kwargs):

    """ Write an sdf file from a dataframe.

    @param df             Pandas object
    @param sdf            A file path provided as a :str:, or a :file-like: object.
    @param write_cols     :bool: specifying whether columns should be written as props
    @param index_as_name  :bool: specifying whether to use index as the name field
    @param mol_props      :bool: specifying whether to write props on the mol in
                          addition to fields in the frame.
    """

    if isinstance(df, pd.Series):
        df = df.to_frame(name='structure')

    writer = Chem.SDWriter(sdf, *args, **kwargs)

    cols = list(df.columns.drop('structure'))

    if not mol_props:
        df.apply(_drop_props, axis=1)

    if write_cols:
        df.apply(_set_props, cols=cols, axis=1)

    if index_as_name:
        df.apply(_set_name, axis=1)

    df.structure.apply(writer.write)


def to_sdf_series(self, *args, **kwargs):

    """ sdf series """

    return write_sdf(self, write_cols=False, *args, **kwargs)


def to_sdf_df(self, *args, **kwargs):

    """ sdf dataframe """

    return write_sdf(self, *args, **kwargs)

pd.Series.to_sdf = to_sdf_series
pd.DataFrame.to_sdf = to_sdf_df


@classmethod
def from_sdf(_, *args, **kwargs):

    """ Create a DataFrame from an sdf file """

    return read_sdf(*args, **kwargs)

pd.DataFrame.from_sdf = from_sdf
