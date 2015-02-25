import pandas as _pd
from rdkit.Chem import DataStructs as _DataStructs
import numpy as _np
import skchem as _skc

def skchemize(func, columns=None, *args, **kwargs):
    """

    transform an RDKit fingerprinting function to work well with pandas

    >>> from rdkit import Chem
    >>> from skchem import *
    >>> f = skchemize(Chem.RDKFingerprint)
    >>> m = Mol.from_smiles('c1ccccc1')
    >>> f(m)
    0     0
    1     0
    2     0
    3     0
    4     0
    5     0
    6     0
    7     0
    8     0
    9     0
    10    0
    11    0
    12    0
    13    0
    14    0
    ...
    2033    0
    2034    0
    2035    0
    2036    0
    2037    0
    2038    0
    2039    0
    2040    0
    2041    0
    2042    0
    2043    0
    2044    0
    2045    0
    2046    0
    2047    0
    Length: 2048, dtype: int64

    >>> df = skchem.read_sdf('skchem/tests/test_resources/hydrocarbons.sdf')
    >>> df.structure.apply(f)
         0     1     2     3     4     5     6     7     8     9     ...   2038  \
    Name                                                              ...          
    297      0     0     0     0     0     0     0     0     0     0  ...      0   
    6324     0     0     0     0     0     0     0     0     0     0  ...      0   
    6334     0     0     0     0     0     0     0     0     0     0  ...      0   

          2039  2040  2041  2042  2043  2044  2045  2046  2047  
    Name                                                        
    297      0     0     0     0     0     0     0     0     0  
    6324     0     0     0     0     0     0     0     0     0  
    6334     0     0     0     0     0     0     0     0     0  

    [3 rows x 2048 columns]
    

    """

    def func_wrapper(m):
        a = _np.array(0)
        _DataStructs.ConvertToNumpyArray(func(m, *args, **kwargs), a)
        return _pd.Series(a, index=columns)
    return func_wrapper

class Fingerprinter(object):

    @classmethod
    def from_rdkit_func(self, func, columns=None, *args, **kwargs):
        fp = Fingerprinter()
        fp.columns = columns
        fp.func = skchemize(func, columns=columns, *args, **kwargs)
        return fp

    def __call__(self, obj):
        return self.calculate(obj)

    def calculate(self, obj):
        if obj.__class__ is _skc.core.Mol:
            return self._calculate_m(obj)

        elif obj.__class__ in [_pd.DataFrame, _pd.Series]:
            return self._calculate_df(obj)

    def _calculate_m(self, m):
        return self.func(m)

    def _calculate_df(self, df):
        return df.structure.apply(self.func)

