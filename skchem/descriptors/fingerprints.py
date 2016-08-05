#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.fingerprints

Fingerprinting classes and associated functions are defined.
"""

import pandas as pd
from rdkit.Chem import GetDistanceMatrix
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
from ..base import Transformer, Featurizer


class MorganFeaturizer(Transformer, Featurizer):
    """ Morgan fingerprints, implemented by RDKit.

     Notes:

         Currently, folded bits are by far the fastest implementation.

     Examples:

        >>> import skchem
        >>> import pandas as pd
        >>> pd.options.display.max_rows = pd.options.display.max_columns = 5

        >>> mf = skchem.descriptors.MorganFeaturizer()
        >>> m = skchem.Mol.from_smiles('CCC')

        Can transform an individual molecule to yield a Series:

        >>> mf.transform(m)
        morgan_fp_idx
        0       0
        1       0
               ..
        2046    0
        2047    0
        Name: MorganFeaturizer, dtype: uint8

        Can transform a list of molecules to yield a DataFrame:

        >>> mf.transform([m])
        morgan_fp_idx  0     1     ...   2046  2047
        0                 0     0  ...      0     0
        <BLANKLINE>
        [1 rows x 2048 columns]

        Change the number of features the fingerprint is folded down to using `n_feats`.

        >>> mf.n_feats = 1024
        >>> mf.transform(m)
        morgan_fp_idx
        0       0
        1       0
               ..
        1022    0
        1023    0
        Name: MorganFeaturizer, dtype: uint8

        Count fingerprints with `as_bits` = False

        >>> mf.as_bits = False
        >>> res = mf.transform(m); res[res > 0]
        morgan_fp_idx
        33     2
        80     1
        294    2
        320    1
        Name: MorganFeaturizer, dtype: int64

        Pseudo-gradient with `grad` shows which atoms contributed to which feature.

        >>> mf.grad(m)[res > 0]
        atom_idx  0  1  2
        features
        33        1  0  1
        80        0  1  0
        294       1  2  1
        320       1  1  1

    """
    def __init__(self, radius=2, n_feats=2048, as_bits=True, use_features=False,
                 use_bond_types=True, use_chirality=False, **kwargs):

        """ Initialize the fingerprinter object.

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

        super(MorganFeaturizer, self).__init__(**kwargs)
        self.radius = radius
        self.n_feats = n_feats
        self.sparse = self.n_feats < 0
        self.as_bits = as_bits
        self.use_features = use_features
        self.use_bond_types = use_bond_types
        self.use_chirality = use_chirality

    def _transform_mol(self, mol):

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
                                               nBits=self.n_feats,
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

    @property
    def name(self):
        return 'morg'

    @property
    def columns(self):
        return pd.RangeIndex(self.n_feats, name='morgan_fp_idx')

    def grad(self, mol):

        """ Calculate the pseudo gradient with respect to the atoms.

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

        cols = pd.Index(list(range(len(mol.atoms))), name='atom_idx')
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


class AtomPairFeaturizer(Transformer, Featurizer):

    """ Atom Pair Fingerprints, implemented by RDKit. """

    def __init__(self, min_length=1, max_length=30, n_feats=2048, as_bits=False,
                 use_chirality=False, **kwargs):

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

        super(AtomPairFeaturizer, self).__init__(**kwargs)
        self.min_length = min_length
        self.max_length = max_length
        self.n_feats = n_feats
        self.sparse = self.n_feats < 0
        self.as_bits = as_bits
        self.use_chirality = use_chirality

    def _transform_mol(self, mol):

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

    @property
    def name(self):
        return 'atom_pair'

    @property
    def columns(self):
        return pd.RangeIndex(self.n_feats, name='ap_fp_idx')


class TopologicalTorsionFeaturizer(Transformer, Featurizer):

    """ Topological Torsion fingerprints, implemented by RDKit. """

    def __init__(self, target_size=4, n_feats=2048, as_bits=False,
                 use_chirality=False, **kwargs):

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
        super(TopologicalTorsionFeaturizer, self).__init__(**kwargs)

    def _transform_mol(self, mol):
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

    @property
    def names(self):
        return 'top_tort'
    @property
    def columns(self):
        return pd.RangeIndex(self.n_feats, name='tt_fp_idx')


class MACCSFeaturizer(Transformer, Featurizer):

    """ MACCS Keys Fingerprints """

    def __init__(self, **kwargs):
        super(MACCSFeaturizer, self).__init__(**kwargs)
        self.n_feats = 166

    def _transform_mol(self, mol):
        return np.array(list(GetMACCSKeysFingerprint(mol)))[1:]

    @property
    def name(self):
        return 'maccs'

    @property
    def columns(self):
        return pd.Index(
            ['ISOTOPE', '103 < ATOMIC NO. < 256', 'GROUP IVA,VA,VIA PERIODS 4-6 (Ge...)', 'ACTINIDE',
             'GROUP IIIB,IVB (Sc...)', 'LANTHANIDE', 'GROUP VB,VIB,VIIB (V...)', 'QAAA@1', 'GROUP VIII (Fe...)',
             'GROUP IIA (ALKALINE EARTH)', '4M RING', 'GROUP IB,IIB (Cu...)', 'ON(C)C', 'S-S', 'OC(O)O', 'QAA@1', 'CTC',
             'GROUP IIIA (B...)', '7M RING', 'SI', 'C=C(Q)Q', '3M RING', 'NC(O)O', 'N-O', 'NC(N)N', 'C$=C($A)$A', 'I',
             'QCH2Q', 'P', 'CQ(C)(C)A', 'QX', 'CSN', 'NS', 'CH2=A', 'GROUP IA (ALKALI METAL)', 'S HETEROCYCLE',
             'NC(O)N', 'NC(C)N', 'OS(O)O', 'S-O', 'CTN', 'F', 'QHAQH', 'OTHER', 'C=CN', 'BR', 'SAN', 'OQ(O)O', 'CHARGE',
             'C=C(C)C', 'CSO', 'NN', 'QHAAAQH', 'QHAAQH', 'OSO', 'ON(O)C', 'O HETEROCYCLE', 'QSQ', 'Snot%A%A', 'S=O',
             'AS(A)A', 'A$A!A$A', 'N=O', 'A$A!S', 'C%N', 'CC(C)(C)A', 'QS', 'QHQH (&...)', 'QQH', 'QNQ', 'NO', 'OAAO',
             'S=A', 'CH3ACH3', 'A!N$A', 'C=C(A)A', 'NAN', 'C=N', 'NAAN', 'NAAAN', 'SA(A)A', 'ACH2QH', 'QAAAA@1', 'NH2',
             'CN(C)C', 'CH2QCH2', 'X!A$A', 'S', 'OAAAO', 'QHAACH2A', 'QHAAACH2A', 'OC(N)C', 'QCH3', 'QN', 'NAAO',
             '5M RING', 'NAAAO', 'QAAAAA@1', 'C=C', 'ACH2N', '8M RING', 'QO', 'CL', 'QHACH2A', 'A$A($A)$A', 'QA(Q)Q',
             'XA(A)A', 'CH3AAACH2A', 'ACH2O', 'NCO', 'NACH2A', 'AA(A)(A)A', 'Onot%A%A', 'CH3CH2A', 'CH3ACH2A',
             'CH3AACH2A', 'NAO', 'ACH2CH2A > 1', 'N=A', 'HETEROCYCLIC ATOM > 1 (&...)', 'N HETEROCYCLE', 'AN(A)A',
             'OCO', 'QQ', 'AROMATIC RING > 1', 'A!O!A', 'A$A!O > 1 (&...)', 'ACH2AAACH2A', 'ACH2AACH2A',
             'QQ > 1 (&...)', 'QH > 1', 'OACH2A', 'A$A!N', 'X (HALOGEN)', 'Nnot%A%A', 'O=A > 1', 'HETEROCYCLE',
             'QCH2A > 1 (&...)', 'OH', 'O > 3 (&...)', 'CH3 > 2 (&...)', 'N > 1', 'A$A!O', 'Anot%A%Anot%A',
             '6M RING > 1', 'O > 2', 'ACH2CH2A', 'AQ(A)A', 'CH3 > 1', 'A!A$A!A', 'NH', 'OC(C)C', 'QCH2A', 'C=O',
             'A!CH2!A', 'NA(A)A', 'C-O', 'C-N', 'O > 1', 'CH3', 'N', 'AROMATIC', '6M RING', 'O', 'RING', 'FRAGMENTS'],
            name='maccs_idx')


class ErGFeaturizer(Transformer, Featurizer):

    """ Extended Reduced Graph Fingerprints.

     Implemented in RDKit."""

    def __init__(self, atom_types=0, fuzz_increment=0.3, min_path=1, max_path=15,  **kwargs):

        super(ErGFeaturizer, self).__init__(**kwargs)
        self.atom_types = atom_types
        self.fuzz_increment = fuzz_increment
        self.min_path = min_path
        self.max_path = max_path
        self.n_feats = 315

    def _transform_mol(self, mol):

        return np.array(GetErGFingerprint(mol))

    @property
    def name(self):
        return 'erg'

    @property
    def columns(self):
        return pd.RangeIndex(self.n_feats, name='erg_fp_idx')


class FeatureInvariantsFeaturizer(Transformer, Featurizer):

    """ Feature invariants fingerprints. """

    def __init__(self, **kwargs):

        super(FeatureInvariantsFeaturizer, self).__init__(**kwargs)

    def _transform_mol(self, mol):

        return np.array(GetFeatureInvariants(mol))

    @property
    def name(self):
        return 'feat_inv'

    @property
    def columns(self):
        return None

class ConnectivityInvariantsFeaturizer(Transformer, Featurizer):

    """ Connectivity invariants fingerprints """

    def __init__(self, include_ring_membership=True, **kwargs):
        super(ConnectivityInvariantsFeaturizer, self).__init__(self, **kwargs)
        self.include_ring_membership = include_ring_membership
        raise NotImplementedError # this is a sparse descriptor

    def _transform_mol(self, mol):

        return np.array(GetConnectivityInvariants(mol))

    @property
    def name(self):
        return 'conn_inv'

    @property
    def columns(self):
        return None

class RDKFeaturizer(Transformer, Featurizer):

    """ RDKit fingerprint """

    # TODO: finish docstring

    def __init__(self, min_path=1, max_path=7, n_feats=2048, n_bits_per_hash=2,
                 use_hs=True, target_density=0.0, min_size=128,
                 branched_paths=True, use_bond_types=True, **kwargs):

        """ RDK fingerprints

        Args:
            min_path (int):
                minimum number of bonds to include in the subgraphs.

            max_path (int):
                maximum number of bonds to include in the subgraphs.

            n_feats (int):
                The number of features to which to fold the fingerprint down. For unfolded, use `-1`.

            n_bits_per_hash (int)
                number of bits to set per path.

            use_hs (bool):
                include paths involving Hs in the fingerprint if the molecule has explicit Hs.

            target_density (float):
                fold the fingerprint until this minimum density has been reached.

            min_size (int):
                the minimum size the fingerprint will be folded to when trying to reach tgtDensity.

            branched_paths (bool):
                if set both branched and unbranched paths will be used in the fingerprint.

            use_bond_types (bool):
                if set both bond orders will be used in the path hashes.

        """

        super(RDKFeaturizer, self).__init__(**kwargs)

        self.min_path = min_path
        self.max_path = max_path
        self.n_feats = n_feats
        self.n_bits_per_hash = n_bits_per_hash
        self.use_hs = use_hs
        self.target_density = target_density
        self.min_size = min_size
        self.branched_paths = branched_paths
        self.use_bond_types = use_bond_types

    def _transform_mol(self, mol):

        return np.array(list(RDKFingerprint(mol, minPath=self.min_path,
                                            maxPath=self.max_path,
                                            fpSize=self.n_feats,
                                            nBitsPerHash=self.n_bits_per_hash,
                                            useHs=self.use_hs,
                                            tgtDensity=self.target_density,
                                            minSize=self.min_size,
                                            branchedPaths=self.branched_paths,
                                            useBondOrder=self.use_bond_types)))

    @property
    def name(self):
        return 'rdkit'

    @property
    def columns(self):
        return pd.RangeIndex(self.n_feats, name='rdk_fp_idx')