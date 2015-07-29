#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.io.smiles

Defining input and output operations for smiles files.
"""

import skchem
import pandas as pd
from skchem.utils import Suppressor

def read_smiles(smiles_file, smiles_column=0, name_column=None, delimiter='\t', \
                        title_line=False, force=False, *args, **kwargs):

    """
    Read a smiles file into a pandas dataframe.  The class wraps the pandas read_csv function.

    @param smiles_file      A file path provided as a :str:, or a :file-like: object.
    @param smiles_column    The column index as an :int: in which the smiles strings are provided.
                            Defaults to _zero_.
    @param name_column      The column index as an :int: in which compound names are provided,
                            for use as the index in the dataframe.  Defaults to _None_.
    @param delimiter        The delimiter used, specified as a :str:.
                            Defaults to _<tab>_.
    @param title_line       A :bool: specifying whether a title line is provided,
                            to use as column titles.
    @param force            A :bool: specifying whether poorly parsed molecules should be skipped,
                            or an error thrown.
    Additionally, pandas read_csv arguments may be provided.

    @returns df             A dataframe of type :pandas.core.frame.DataFrame:.
    """

    with Suppressor():

        # set the header line to pass to the pandas parser
        # we accept True as being line zero, as is usual for smiles
        # if user specifies a header already, then do nothing

        header = kwargs.get('header', None)
        if title_line is True:
            header = 0
        elif header is not None:
            kwargs.pop('header') #remove from the kwargs to not pass it twice
        else:
            header = None

        # open file if not already open
        if isinstance(smiles_file, str):
            smiles_file = open(smiles_file, 'r')

        # read the smiles file
        df = pd.read_csv(smiles_file, delimiter=delimiter, header=header, *args, **kwargs)

        # replace the smiles column with the structure column
        lst = list(df.columns)
        lst[smiles_column] = 'structure'
        df.columns = lst

        # apply the from smiles constructor
        if force:
            def parse(smiles):

                """
                Parse a molecule from smiles string and return None if it doesn't load
                (restoring rdkit functionality)
                """

                try:
                    return skchem.Mol.from_smiles(smiles)
                except ValueError:
                    return None
        else:
            def parse(smiles):
                """
                Parse a molecule from smiles string
                """
                return skchem.Mol.from_smiles(smiles)

        df['structure'] = df['structure'].apply(str).apply(parse) #make sure is a string

        if force:
            df = df[df['structure'].notnull()]

        # set index if passed
        if name_column is not None:
            df = df.set_index(df.columns[name_column])

        return df

@classmethod
def from_smiles(_, *args, **kwargs):
    """ Create a DataFrame from a smiles file """
    return read_smiles(*args, **kwargs)

#set on pandas dataframe
pd.DataFrame.from_smiles = from_smiles
