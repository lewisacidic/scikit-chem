#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.fingerprints

Fingerprinting classes and associated functions are defined.
"""

from functools import wraps

import pandas as pd
from rdkit.Chem import DataStructs, GetDistanceMatrix
from rdkit.Chem.rdMolDescriptors import (GetMorganFingerprint,
                                         GetHashedMorganFingerprint,
                                         GetAtomPairFingerprint,
                                         GetHashedAtomPairFingerprint,
                                         GetTopologicalTorsionFingerprint,
                                         GetHashedTopologicalTorsionFingerprint,
                                         GetMACCSKeysFingerprint,
                                         GetFeatureInvariants,
                                         GetConnectivityInvariants)
from rdkit.Chem.rdReducedGraphs import GetErGFingerprint
from rdkit.Chem.rdmolops import RDKFingerprint

import numpy as np
import skchem

def skchemize(func, columns=None, *args, **kwargs):
    """

    transform an RDKit fingerprinting function to work well with pandas

    >>> from rdkit import Chem
    >>> import skchem
    >>> from skchem.descriptors import skchemize
    >>> from skchem.core import Mol
    >>> f = skchemize(Chem.RDKFingerprint)
    >>> m = Mol.from_smiles('c1ccccc1')
    >>> f(m)
    0       0
    1       0
    2       0
    3       0
    4       0
    5       0
    6       0
    7       0
    8       0
    9       0
    10      0
    11      0
    12      0
    13      0
    14      0
    15      0
    16      0
    17      0
    18      0
    19      0
    20      0
    21      0
    22      0
    23      0
    24      0
    25      0
    26      0
    27      0
    28      0
    29      0
           ..
    2018    0
    2019    0
    2020    0
    2021    0
    2022    0
    2023    0
    2024    0
    2025    0
    2026    0
    2027    0
    2028    0
    2029    0
    2030    0
    2031    0
    2032    0
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
    dtype: int64
    >>> from skchem.data import resource
    >>> df = skchem.read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
    >>> df.structure.apply(f)
          0     1     2     3     4     5     6     7     8     9     ...   2038
    name                                                              ...
    297      0     0     0     0     0     0     0     0     0     0  ...      0
    6324     0     0     0     0     0     0     0     0     0     0  ...      0
    6334     0     0     0     0     0     0     0     0     0     0  ...      0
    <BLANKLINE>
          2039  2040  2041  2042  2043  2044  2045  2046  2047
    name
    297      0     0     0     0     0     0     0     0     0
    6324     0     0     0     0     0     0     0     0     0
    6334     0     0     0     0     0     0     0     0     0
    <BLANKLINE>
    [3 rows x 2048 columns]

    """
    @wraps(func)
    def func_wrapper(m):

        """ Function that wraps an rdkit function allowing it to produce dataframes. """

        arr = np.array(0)
        DataStructs.ConvertToNumpyArray(func(m, *args, **kwargs), arr)

        return pd.Series(arr, index=columns)

    return func_wrapper


class Fingerprinter(object):

    """ Fingerprinter Base class. """

    def __init__(self, func, name=None):

        """ A generic fingerprinter.  Create with a function.

        Args:
            func (callable):
                A fingerprinting function that takes an skchem.Mol argument, and
                returns an iterable of values.
            name (str):
                The name of the fingerprints that are being calculated"""

        self.NAME = name
        self.func = func

    def __call__(self, obj):

        """ Call the fingerprinter directly.

        This is a shorthand for transform. """

        return self.transform(obj)

    def __add__(self, other):

        """ Add fingerprinters together to create a fusion fingerprinter.

        Fusion featurizers will transform molecules to series with all
        features from all component featurizers.
        """

        fpers = []
        for fper in (self, other):
            if isinstance(fper, FusionFingerprinter):
                fpers += fper.fingerprinters
            else:
                fpers.append(fper)

        return FusionFingerprinter(fpers)

    def fit(self, X, y):
        return self

    def transform(self, obj):

        """ calculate the fingerprint for the given object. """

        if isinstance(obj, skchem.Mol):
            return self._transform(obj)

        elif isinstance(obj, pd.DataFrame):
            return obj.structure.apply(self.transform)

        elif isinstance(obj, pd.Series):
            return obj.apply(self.transform)

        else:
            raise NotImplementedError

    def _transform(self, mol):

        """ Calculate the fingerprint on a molecule. """

        return pd.Series(list(self.func(mol)), name=mol.name)


class FusionFingerprinter(Fingerprinter):

    def __init__(self, fingerprinters):

        self.fingerprinters = fingerprinters

    def transform(self, obj):

        if isinstance(obj, skchem.Mol):
            return pd.concat([fp.transform(obj) for fp in self.fingerprinters],
                             keys=[fp.NAME for fp in self.fingerprinters])

        elif isinstance(obj, pd.DataFrame):
            return pd.concat([fp.transform(obj) for fp in self.fingerprinters],
                             keys=[fp.NAME for fp in self.fingerprinters],
                             axis=1)

        elif isinstance(obj, pd.Series):
            return pd.concat([fp.transform(obj.structure) \
                                for fp in self.fingerprinters],
                             keys=[fp.NAME for fp in self.fingerprinters],
                             axis=1)

        else:
            raise NotImplementedError

    def _transform(self, mol):

        return pd.concat([fp.transform(mol) for fp in self.fingerprinters])

class MorganFingerprinter(Fingerprinter):

    """ Morgan Fingerprint Transformer. """

    NAME = 'morgan'

    def __init__(self, radius=2, n_feats=2048, as_bits=False,
                 use_features=False, use_bond_types=True, use_chirality=False):

        """
        Args:
            radius (int):
                The maximum radius for atom environments.
                Default is `2`.
            n_feats (int):
                The number of features to which to fold the fingerprint down.
                For unfolded, use `-1`.
                Default is `2048`.
            as_bits (bool):
                Whether to return bits (`True`) or counts (`False`).
                Default is `True`.
            use_features (bool):
                Whether to use map atom types to generic features (FCFP analog).
                Default is `False`.
            use_bond_types (bool):
                Whether to use bond types to differentiate environments.
                Default is `False`.
            use_chirality (bool):
                Whether to use chirality to differentiate environments.
                Default is `False`.
        """

        self.radius = radius
        self.n_feats = n_feats
        self.as_bits = as_bits
        self.use_features = use_features
        self.use_bond_types = use_bond_types
        self.use_chirality = use_chirality

    def _transform(self, mol):

        """Private method to transform a skchem molecule.

        Use transform` for the public method, which genericizes the argument to
        iterables of mols.

        Args:
            mol (skchem.Mol): Molecule to calculate fingerprint for.

        Returns:
            pandas.Series:
                Fingerprint as a series.
        """

        if self.n_feats < 0:

            res = GetMorganFingerprint(mol, self.radius,
                                       useFeatures=self.use_features,
                                       useBondTypes=self.use_bond_types,
                                       useChirality=self.use_chirality).GetNonzeroElements()
            idx = pd.Index(list(res.keys()), name='features')

        else:
            res = list(GetHashedMorganFingerprint(mol, self.radius,
                                        nBits=self.n_feats,
                                        useFeatures=self.use_features,
                                        useBondTypes=self.use_bond_types,
                                        useChirality=self.use_chirality))
            idx = pd.Index(range(self.n_feats), name='features')

        res = pd.Series(res, name=mol.name, index=idx).astype(int)

        if self.as_bits:
            res = (res > 0).astype(np.uint8) #smallest memory size

        return res

    def grad(self, mol):

        """ Calculate the pseudo gradient with resepect to the atoms.

        The pseudo gradient is the number of times the atom set that particular
        bit.

        Args:
            mol (skchem.Mol):
                The molecule for which to calculate the pseudo gradient.

        Returns:
            pandas.DataFrame:
                Dataframe of pseudogradients, with columns corresponding to
                atoms, and rows corresponding to features of the fingerprint.
        """

        cols = pd.Index(list(range(len(mol.atoms))), name='atoms')
        dist = GetDistanceMatrix(mol)

        info = {}

        if self.n_feats < 0:

            res = GetMorganFingerprint(mol, self.radius,
                                       useFeatures=self.use_features,
                                       useBondTypes=self.use_bond_types,
                                       useChirality=self.use_chirality,
                                       bitInfo=info).GetNonzeroElements()
            idx_list = list(res.keys())
            idx = pd.Index(idx_list, name='features')
            grad = np.zeros((len(idx), len(cols)))
            for bit in info:
                for atom_idx, radius in info[bit]:
                    grad[idx_list.index(bit)] += (dist <= radius)[atom_idx]

        else:

            res = list(GetHashedMorganFingerprint(mol, self.radius,
                                        nBits=self.n_feats,
                                        useFeatures=self.use_features,
                                        useBondTypes=self.use_bond_types,
                                        useChirality=self.use_chirality,
                                        bitInfo=info))
            idx = pd.Index(range(self.n_feats), name='features')
            grad = np.zeros((len(idx), len(cols)))

            for bit in info:
                for atom_idx, radius in info[bit]:
                    grad[bit] += (dist <= radius)[atom_idx]

        grad = pd.DataFrame(grad, index=idx, columns=cols)

        if self.as_bits:
            grad = (grad > 0)

        return grad.astype(int)

class AtomPairFingerprinter(Fingerprinter):

    """ Atom Pair Tranformer. """

    NAME = 'atom_pair'

    def __init__(self, min_length=1, max_length=30, n_feats=2048, as_bits=False,
                 use_chirality=False):

        """ Instantiate an atom pair fingerprinter.

        Args:
            min_length (int):
                The minimum length of paths between pairs.
                Default is `1`, i.e. pairs can be bonded together.
            max_length (int):
                The maximum length of paths between pairs.
                Default is `30`.
            n_feats (int):
                The number of features to which to fold the fingerprint down.
                For unfolded, use `-1`.
                Default is `2048`.
            as_bits (bool):
                Whether to return bits (`True`) or counts (`False`).
                Default is `False`.
            use_chirality (bool):
                Whether to use chirality to differentiate environments.
                Default is `False`.
        """

        self.min_length = min_length
        self.max_length = max_length
        self.n_feats = n_feats
        self.as_bits = as_bits
        self.use_chirality = use_chirality

    def _transform(self, mol):

        """Private method to transform a skchem molecule.

        Use transform` for the public method, which genericizes the argument to
        iterables of mols.

        Args:
            mol (skchem.Mol): Molecule to calculate fingerprint for.

        Returns:
            pandas.Series:
                Fingerprint as a series.
        """

        if self.n_feats == -1:

            res = GetAtomPairFingerprint(mol, minLength=self.min_length,
                                         maxLength=self.max_length,
                                         includeChirality=self.use_chirality)
        else:
            res = list(GetHashedAtomPairFingerprint(mol,
                                        minLength=self.min_length,
                                         maxLength=self.max_length,
                                         nBits=self.n_feats,
                                         includeChirality=self.use_chirality))

        res = pd.Series(res, name=mol.name)

        if self.as_bits:
            return (res > 0).astype(int)
        else:
            return res

class TopologicalTorsionFingerprinter(Fingerprinter):

    NAME = 'topological_torsion'

    def __init__(self, target_size=4, n_feats=2048, as_bits=False,
                 use_chirality=False):

        """
        Args:
            target_size (int):
                # TODO
            n_feats (int):
                The number of features to which to fold the fingerprint down.
                For unfolded, use `-1`.
                Default is `2048`.
            as_bits (bool):
                Whether to return bits (`True`) or counts (`False`).
                Default is `False`.
            use_chirality (bool):
                Whether to use chirality to differentiate environments.
                Default is `False`.
        """

        self.target_size = target_size
        self.n_feats = n_feats
        self.as_bits = as_bits
        self.use_chirality = use_chirality

    def _transform(self, mol):

        if self.n_feats == -1:

            res = GetTopologicalTorsionFingerprint(mol,
                                        targetSize=self.targetSize,
                                        includeChirality=self.use_chirality)

        else:
            res = list(GetHashedTopologicalTorsionFingerprint(mol,
                                        targetSize=self.targetSize,
                                        nBits=self.n_feats))

        res = pd.Series(res, name=mol.name)

        if self.as_bits:
            return (res > 0).astype(int)
        else:
            return res


class MACCSKeysFingerprinter(Fingerprinter):

    """ MACCS Keys Fingerprints """

    NAME = 'maccs'

    def __init__(self):
        pass

    def _transform(self, mol):

        return pd.Series(list(GetMACCSKeysFingerprint(mol)))

class ErGFingerprinter(Fingerprinter):

    """ ErG Fingerprints """

    NAME = 'erg'

    def __init__(self):
        pass

    def _transform(self, mol):

        return pd.Series(GetErGFingerprint(mol))

class FeatureInvariantsFingerprinter(Fingerprinter):

    """ Feature invariants fingerprints. """

    NAME = 'feat_inv'

    def __init__(self):
        pass

    def _transform(self, mol):

        return pd.Series(GetFeatureInvariants(mol))

class ConnectivityInvariantsFingerprinter(Fingerprinter):

    """ Connectivity invariants fingerprints """

    NAME = 'conn_inv'

    def __init__(self):
        pass

    def _transform(self, mol):

        return pd.Series(GetConnectivityInvariants(mol))

class RDKFingerprinter(Fingerprinter):

    """ RDKit fingerprint """

    NAME = 'rdk'

    def __init__(self, min_path=1, max_path=7, n_feats=2048, n_bits_per_hash=2,
                 use_hs=True, target_density=0.0, min_size=128,
                 branched_paths=True, use_bond_types=True):

        """ RDK fingerprints

        # TODO

        Args:
            min_path (int):

            max_path (int):

            n_feats (int):
                The number of features to which to fold the fingerprint down.
                For unfolded, use `-1`.
                Default is `2048`.

            n_bits_per_hash (int)

            use_hs (bool):

            target_density (float):

            min_size (int):

            branched_paths (bool):

            use_bond_types (bool):
        """

        self.min_path = 1
        self.max_path = 7
        self.n_feats = 2048
        self.n_bits_per_hash = 2
        self.use_hs = True
        self.target_density = 0.0
        self.min_size = 128
        self.branched_paths = True
        self.use_bond_types = True

    def _transform(self, mol):

        return pd.Series(list(RDKFingerprint(mol, minPath=self.min_path,
                                             maxPath=self.max_path,
                                             fpSize=self.n_feats,
                                             nBitsPerHash=self.n_bits_per_hash,
                                             useHs=self.use_hs,
                                             tgtDensity=self.target_density,
                                             minSize=self.min_size,
                                             branchedPaths=self.branched_paths,
                                             useBondOrder=self.use_bond_types)),
                        name=mol.name)
