#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""skchem.io.smiles

Defining input and output operations for smiles files."""

from rdkit import Chem
import skchem
import pandas as pd

def read_smiles(smiles_file, smiles_column=0, name_column=None, delimiter='\t', title_line=True, *args, **kwargs):

    """
    Read a smiles file into a pandas dataframe.  The class wraps the pandas read_csv function.

    @param smiles_file      A file path provided as a :str:, or a :file-like: object.
    @param smiles_column    The column index as an :int: in which the smiles strings are provided. 
                            Defaults to _zero_.
    @param name_column      The column index as an :int: in which compound names are provided, for use as the index in the dataframe.
                            Defaults to _None_.
    @param delimiter        The delimiter used, specified as a :str:.  
                            Defaults to _<tab>_.
    @param title_line       A :bool: specifying whether a title line is provided, to use as column titles.

    Additionally, pandas read_csv arguments may be provided.

    @returns df             A dataframe of type :pandas.core.frame.DataFrame:.
    """

    # set the header line to pass to the pandas parser
    # we accept True as being line zero, as is usual for smiles
    if title_line is True:
        header = 0
    else:
        header = None

    # open file if not already open
    if type(smiles_file) is str:
        smiles_file = open(smiles_file, 'r')

    # read the smiles file
    df = pd.read_csv(smiles_file, delimiter=delimiter, header=header)

    # replace the smiles column with the structure column
    l = list(df.columns)
    l[smiles_column] = 'structure'
    df.columns = l

    # apply the from smiles constructor
    df['structure'] = df['structure'].apply(lambda smi: skchem.Mol.from_smiles(smi)) 

    # set index if passed
    if name_column is not None:
        df = df.set_index(df.columns[name_column])
    
    return df

@classmethod
def from_smiles(self, *args, **kwargs):
    return read_smiles(*args, **kwargs)

#set on pandas dataframe
pd.DataFrame.from_smiles = from_smiles