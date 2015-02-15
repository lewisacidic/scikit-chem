from rdkit import Chem as _Chem
import pandas as _pd

def read_sdf(sdf_file, *args, **kwargs):

    if type(sdf_file) is str:
        sdf_file = open(sdf_file, 'r')

    ms = []
    idx = []
    props = set()

    mol_supp = _Chem.ForwardSDMolSupplier(sdf_file, *args, **kwargs)
    for m in mol_supp:
        if m is None:
            continue
        ms.append(m)
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

    df.index.name = 'Name'

    return df

def read_smiles(smiles_file, *args, **kwargs):

    if type(smiles_file) is str:
        smiles_file = open(smiles_file, 'r')

    df = _pd.read_csv(smiles_file, delimiter='\t', header=0)

    df.columns = [col if i > 0 else 'structure' for i, col in enumerate(df.columns)]
    df['structure'] = df['structure'].apply(_Chem.MolFromSmiles)

    return df

@classmethod
def from_sdf(self, *args, **kwargs):
    return read_sdf(*args, **kwargs)

#set on pandas dataframe
_pd.DataFrame.from_sdf = from_sdf