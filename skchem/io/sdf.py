#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.io.sdf

Defining input and output operations for sdf files.
"""

from rdkit import Chem
import skchem
from skchem.utils import Suppressor
import pandas as pd

def read_sdf(sdf_file, force=False, *args, **kwargs):

    """
        Read an sdf file into a pandas dataframe.
        The function wraps the RDKit ForwardSDMolSupplier function.

        @param sdf_file     A file path provided as a :str:, or a :file-like: object.
        @param force        A :bool: specifying if the dataframe should be constructed even if
                            a molecule fails to parse correctly.  Defaults to _True_.

        Additionally, ForwardSDMolSupplier arguments may be passed.

        @returns df         A dataframe of type :pandas.core.frame.DataFrame:.

    """

    # use the suppression context manager to not pollute our stdout with rdkit errors and warnings.
    # perhaps this should be captured better by Mol etc.

    with Suppressor():

        if isinstance(sdf_file, str):
            sdf_file = open(sdf_file, 'rb') # use read bytes for python 3 compatibility

        # use a list as we need a resizing array as we don't know how many compounds there are.
        ms = []
        idx = []
        props = set()

        mol_supp = Chem.ForwardSDMolSupplier(sdf_file, *args, **kwargs)
        for i, m in enumerate(mol_supp):
            if m is None and force is False:
                # raise Value Error, like in the json module when no json is detected
                raise ValueError('Molecule {} could not be decoded.'.format(i + 1))
            elif m is None and force is True:
                continue
            ms.append(skchem.Mol(m))
            idx.append(m.GetProp('_Name'))
            props = props.union(props, set(m.GetPropNames()))

        df = pd.DataFrame(ms, index=idx, columns=['structure'])

        for prop in props:
            df[prop] = df.structure.apply(lambda m: m.props.get(prop, None))

        df.index.name = 'name'

        return df

@classmethod
def from_sdf(_, *args, **kwargs):

    """ Create a DataFrame from an sdf file """

    return read_sdf(*args, **kwargs)

#set on pandas dataframe
pd.DataFrame.from_sdf = from_sdf
