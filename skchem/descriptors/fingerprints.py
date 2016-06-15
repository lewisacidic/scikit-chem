#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.fingerprints

Fingerprinting classes and associated functions are defined.
"""

from functools import wraps
from collections import Iterable
import pandas as pd
from rdkit.Chem import DataStructs, GetDistanceMatrix
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit.Chem.rdMolDescriptors import (GetMorganFingerprint,
                                         GetHashedMorganFingerprint,
                                         GetMorganFingerprintAsBitVect,
                                         GetAtomPairFingerprint,
                                         GetHashedAtomPairFingerprint,
                                         GetHashedAtomPairFingerprintAsBitVect,
                                         GetTopologicalTorsionFingerprint,
                                         GetHashedTopologicalTorsionFingerprint,
                                         GetHashedTopologicalTorsionFingerprintAsBitVect,
                                         GetMACCSKeysFingerprint,
                                         GetFeatureInvariants,
                                         GetConnectivityInvariants)
from rdkit.Chem.rdReducedGraphs import GetErGFingerprint
from rdkit.Chem.rdmolops import RDKFingerprint

import numpy as np
import skchem

class Fingerprinter(object):

    """ Fingerprinter Base class. """

    def __init__(self, func, sparse=False, name=None):

        """ A generic fingerprinter.  Create with a function.

        Args:
            func (callable):
                A fingerprinting function that takes an skchem.Mol argument, and
                returns an iterable of values.
            name (str):
                The name of the fingerprints that are being calculated"""

        self.NAME = name
        self.func = func
        self.sparse = sparse

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

        """ Calculate a fingerprint for the given object.

        Args:
            obj (skchem.Mol or pd.Series or pd.DataFrame or iterable):
                The object to be transformed.

        Returns:
            pd.DataFrame:
                The produced features.
        """

        if self.sparse:
            return self._transform_sparse(obj)
        else:
            return self._transform_dense(obj)

    def _transform_dense(self, obj):

        """ calculate a dense fingerprint for the given object. """

        if isinstance(obj, skchem.Mol):
            return pd.Series(self._transform(obj), index=self.index)

        elif isinstance(obj, pd.DataFrame):
            return self.transform(obj.structure)

        elif isinstance(obj, pd.Series):
            res_0 = self._transform(obj.iloc[0])
            res = np.zeros((len(obj), len(res_0)))
            for i, mol in enumerate(obj):
                res[i] = self._transform(mol)
            return pd.DataFrame(res, index=obj.index, columns=self.index)

        elif isinstance(obj, (tuple, list)):
            res_0 = self._transform(obj[0])
            res = np.zeros((len(obj), len(res_0)))
            for i, mol in enumerate(obj):
                res[i] = self._transform(mol)

            idx = pd.Index([mol.name for mol in obj], name='name')
            return pd.DataFrame(res, index=idx, columns=self.index)

        else:
            raise NotImplementedError

    def _transform_sparse(self, obj):

        """ Calculate a sparse fingerprint for the given object. """

        if isinstance(obj, skchem.Mol):
            return pd.Series(self._transform(obj), index=self.index)

        elif isinstance(obj, pd.DataFrame):
            return self.transform(obj.structure)

        elif isinstance(obj, pd.Series):
            return pd.DataFrame([self.transform(m) for m in obj],
                                index=obj.idx,
                                columns=self.index).fillna(0)

        elif isinstance(obj, (tuple, list)):
            idx = pd.Index([mol.name for mol in obj], name='name')
            return pd.DataFrame([self.transform(m) for m in obj],
                                index=idx,
                                columns=self.index)

        else:
            raise NotImplementedError

    def _transform(self, mol):

        """ Calculate the fingerprint on a molecule. """

        return pd.Series(list(self.func(mol)), name=mol.name)

    @property
    def index(self):

        """ The index to use. """

        return None

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

    def __init__(self, radius=2, n_feats=2048, as_bits=True,
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

        Notes:
            Currently, folded bits are by far the fastest implementation.
        """

        self.radius = radius
        self.n_feats = n_feats
        self.sparse = self.n_feats < 0
        self.as_bits = as_bits
        self.use_features = use_features
        self.use_bond_types = use_bond_types
        self.use_chirality = use_chirality

    def _transform(self, mol):

        """Private method to transform a skchem molecule.

        Use `transform` for the public method, which genericizes the argument to
        iterables of mols.

        Args:
            mol (skchem.Mol): Molecule to calculate fingerprint for.

        Returns:
            np.array or dict:
                Fingerprint as an array (or a dict if sparse).
        """

        if self.as_bits and self.n_feats > 0:

            fp = GetMorganFingerprintAsBitVect(mol, self.radius,
                                           useFeatures=self.use_features,
                                           useBondTypes=self.use_bond_types,
                                           useChirality=self.use_chirality)
            res = np.array(0)
            ConvertToNumpyArray(fp, res)
            res = res.astype(np.uint8)

        else:

            if self.n_feats <= 0:

                res = GetMorganFingerprint(mol, self.radius,
                                           useFeatures=self.use_features,
                                           useBondTypes=self.use_bond_types,
                                           useChirality=self.use_chirality)
                res = res.GetNonzeroElements()
                if self.as_bits:
                    res = {k: int(v > 0) for k, v in res.items()}

            else:
                res = GetHashedMorganFingerprint(mol, self.radius,
                                                 nBits=self.n_feats,
                                                 useFeatures=self.use_features,
                                                 useBondTypes=self.use_bond_types,
                                                 useChirality=self.use_chirality)
                res = np.array(list(res))



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
        self.sparse = self.n_feats < 0
        self.as_bits = as_bits
        self.use_chirality = use_chirality

    def _transform(self, mol):

        """Private method to transform a skchem molecule.

        Use transform` for the public method, which genericizes the argument to
        iterables of mols.

        Args:
            mol (skchem.Mol): Molecule to calculate fingerprint for.

        Returns:
            np.array or dict:
                Fingerprint as an array (or a dict if sparse).
        """


        if self.as_bits and self.n_feats > 0:

            fp = GetHashedAtomPairFingerprintAsBitVect(mol, nBits=self.n_feats,
                                           minLength=self.min_length,
                                           maxLength=self.max_length,
                                           includeChirality=self.use_chirality)
            res = np.array(0)
            ConvertToNumpyArray(fp, res)
            res = res.astype(np.uint8)

        else:

            if self.n_feats <= 0:

                res = GetAtomPairFingerprint(mol, nBits=self.n_feats,
                                               minLength=self.min_length,
                                               maxLength=self.max_length,
                                               includeChirality=self.use_chirality)
                res = res.GetNonzeroElements()
                if self.as_bits:
                    res = {k: int(v > 0) for k, v in res.items()}

            else:
                res = GetHashedAtomPairFingerprint(mol, nBits=self.n_feats,
                                               minLength=self.min_length,
                                               maxLength=self.max_length,
                                               includeChirality=self.use_chirality)
                res = np.array(list(res))

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
        self.sparse = self.n_feats < 0
        self.as_bits = as_bits
        self.use_chirality = use_chirality

    def _transform(self, mol):
        """ Private method to transform a skchem molecule.
        Args:
            mol (skchem.Mol): Molecule to calculate fingerprint for.

        Returns:
            np.array or dict:
                Fingerprint as an array (or a dict if sparse).
        """

        if self.as_bits and self.n_feats > 0:

            fp = GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=self.n_feats,
                                           targetSize=self.target_size,
                                           includeChirality=self.use_chirality)
            res = np.array(0)
            ConvertToNumpyArray(fp, res)
            res = res.astype(np.uint8)

        else:

            if self.n_feats <= 0:

                res = GetTopologicalTorsionFingerprint(mol, nBits=self.n_feats,
                                               targetSize=self.target_size,
                                               includeChirality=self.use_chirality)
                res = res.GetNonzeroElements()
                if self.as_bits:
                    res = {k: int(v > 0) for k, v in res.items()}

            else:
                res = GetHashedTopologicalTorsionFingerprint(mol, nBits=self.n_feats,
                                               targetSize=self.target_size,
                                               includeChirality=self.use_chirality)
                res = np.array(list(res))

        return res


class MACCSKeysFingerprinter(Fingerprinter):

    """ MACCS Keys Fingerprints """

    NAME = 'maccs'

    def __init__(self):
        self.sparse = False

    def _transform(self, mol):

        return np.array(list(GetMACCSKeysFingerprint(mol)))

class ErGFingerprinter(Fingerprinter):

    """ ErG Fingerprints """

    NAME = 'erg'

    def __init__(self):
        self.sparse = False

    def _transform(self, mol):

        return np.array(GetErGFingerprint(mol))

class FeatureInvariantsFingerprinter(Fingerprinter):

    """ Feature invariants fingerprints. """

    NAME = 'feat_inv'

    def __init__(self):
        self.sparse = False

    def _transform(self, mol):

        return np.array(GetFeatureInvariants(mol))

class ConnectivityInvariantsFingerprinter(Fingerprinter):

    """ Connectivity invariants fingerprints """

    NAME = 'conn_inv'

    def __init__(self):
        self.sparse = False

    def _transform(self, mol):

        return np.array(GetConnectivityInvariants(mol))

class RDKFingerprinter(Fingerprinter):

    """ RDKit fingerprint """

    NAME = 'rdk'

    def __init__(self, min_path=1, max_path=7, n_feats=2048, n_bits_per_hash=2,
                 use_hs=True, target_density=0.0, min_size=128,
                 branched_paths=True, use_bond_types=True):

        """ RDK fingerprints

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
        self.sparse = False
        self.n_bits_per_hash = 2
        self.use_hs = True
        self.target_density = 0.0
        self.min_size = 128
        self.branched_paths = True
        self.use_bond_types = True

    def _transform(self, mol):

        return np.array(list(RDKFingerprint(mol, minPath=self.min_path,
                                             maxPath=self.max_path,
                                             fpSize=self.n_feats,
                                             nBitsPerHash=self.n_bits_per_hash,
                                             useHs=self.use_hs,
                                             tgtDensity=self.target_density,
                                             minSize=self.min_size,
                                             branchedPaths=self.branched_paths,
                                             useBondOrder=self.use_bond_types)),
                        name=mol.name)
