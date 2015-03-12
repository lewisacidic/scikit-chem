from rdkit import Chem as _Chem
import skchem as _skc
import pandas as _pd
import json
from collections import defaultdict

def read_sdf(sdf_file, *args, **kwargs):

    if type(sdf_file) is str:
        sdf_file = open(sdf_file, 'rb')

    ms = []
    idx = []
    props = set()

    mol_supp = _Chem.ForwardSDMolSupplier(sdf_file, *args, **kwargs)
    for i, m in enumerate(mol_supp):
        if m is None:
            # raise Value Error, like in the json module when no json is detected
            raise ValueError('Molecule {} could not be decoded.'.format(i + 1))
        ms.append(_skc.Mol(m))
        idx.append(m.GetProp('_Name'))
        props = props.union(props, set(m.GetPropNames()))

    df = _pd.DataFrame(ms, idx, ['structure'])

    def get_prop(m, prop):
        '''get the properties for a molecule'''
        try:
            return m.GetProp(prop)
        except KeyError:
            return None

    for prop in props:
        df[prop] = df.structure.apply(lambda m: get_prop(m, prop))

    df.index.name = 'name'

    return df

def read_smiles(smiles_file, smiles_column=0, index_column=None, delimiter='\t', header=True, *args, **kwargs):

    # set the header line to pass to the pandas parser

    if header:
        header = 0
    else:
        header = None

    # open file if not already open
    if type(smiles_file) is str:
        smiles_file = open(smiles_file, 'r')

    # read the smiles file as a ?del
    df = _pd.read_csv(smiles_file, delimiter=delimiter, header=header)

    # replace the smiles column with the structure column 
    l = list(df.columns)
    l[smiles_column] = 'structure'
    df.columns = l

    # apply the from smiles constructor
    df['structure'] = df['structure'].apply(lambda smi: _skc.Mol.from_smiles(smi)) #perhaps should raise errors for incorrectly parsed molecules here?
    
    # set index if passed
    if index_column is not None:
        df = df.set_index(df.columns[index_column])
    
    return df

_to_dict = _pd.DataFrame.to_dict
_to_json = _pd.DataFrame.to_json

def to_dict(self, *args, **kwargs):
    kwargs = defaultdict(lambda: None, kwargs)
    if kwargs['orient'] == 'chemdoodle':
        #this may be wrong - may have to extract
        return {'m': [m.to_dict(kind='chemdoodle')['m'][0] for m in self.structure]}
    else:
        return _to_dict(self, *args, **kwargs)

def to_json(self, *args, **kwargs):
    kwargs = defaultdict(lambda: None, kwargs)
    if kwargs['orient'] == 'chemdoodle':
        return json.dumps(self.to_dict(*args, **kwargs))
    else:
        return _to_dict(self, *args, **kwargs)

_pd.DataFrame.to_dict = to_dict
_pd.DataFrame.to_json = to_json

@classmethod
def _from_sdf(self, *args, **kwargs):
    return read_sdf(*args, **kwargs)

@classmethod
def _from_smiles(self, *args, **kwargs):
    return read_smiles(*args, **kwargs)

#set on pandas dataframe
_pd.DataFrame.from_sdf = _from_sdf
_pd.DataFrame.from_smiles = _from_smiles