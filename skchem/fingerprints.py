import pandas as _pd

def skchemize(func, *args, **kwargs):
    '''transform a fingerprinting function to work well with pandas

    >>> from rdkit import Chem
    >>> import skchem
    >>> f = skchemize(Chem.RDKFingerprint)
    >>> m = Chem.MolFromSmiles('c1ccccc1')
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

    >>> @skchemize
    >>> def weird_fingerprint(m):
    ...     return [m.GetNumAtoms(), m.GetNumBonds()]
    ... 
    >>> df.structure.apply(weird_fingerprint)
          0  1
    Name      
    297   1  0
    6324  2  1
    6334  3  2

    '''

    def func_wrapper(m):
        return _pd.Series(list(func(m, *args, **kwargs)))
    return func_wrapper

