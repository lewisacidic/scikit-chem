#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.cross_validation.similarity_threshold

Similarity threshold dataset partitioning functionality.
"""


import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform, cdist
from scipy.sparse import dok_matrix, triu

from .. import descriptors


class SimThresholdSplit(object):

    def __init__(self, inp, threshold=0.5, fper='morgan',
                 similarity_metric='jaccard', memory_optimized=False,
                 fingerprints=None, similarity_matrix=None):
        """ Threshold similarity split for chemical datasets.

        This class implements a splitting technique that will pool compounds
        with similarity above a theshold into the same splits.

        Machine learning techniques should be able to extrapolate outside of a
        molecular series, or scaffold, however random splits will result in some
        'easy' test sets that are either *identical* or in the same molecular
        series or share a significant scaffold with training set compounds.

        This splitting technique reduces or eliminates (depending on the
        threshold set) this effect, making the problem harder.

        Args:
            inp (scipy.sparse.dok, pd.Series or pd.DataFrame):
                Either:
                - a series of skchem.Mols
                - dataframe of precalculated fingerprints

            n_splits (int):
                The number of splits to give.  This will be overridden if ratio
                is passed.

            ratio (list[floats]):
                Split ratios to use.

            threshold (float):
                The similarity threshold, above which, compounds will all be
                assigned to the same split.

            fper (str or skchem.Fingerprinter):
                The fingerprinting technique to use to generate the similarity
                matrix.

            fingerprints (bool):
                Whether percalculated fingerprints were passed directly.

            similarity_matrix (scipy.sparse.dok):
                A precalculated similarity matrix.

        Notes:
            The splits will not always be exactly the size requested, due to the
            constraint and requirement to maintain random shuffling.
        """

        if isinstance(fper, str):
            fper = descriptors.get(fper)

        self.fper = fper
        fps = inp if fingerprints else self.fper.transform(inp)

        self.n_instances = len(inp)

        self.threshold = threshold
        self.similarity_metric = similarity_metric
        self.memory_optimized = memory_optimized

        if not similarity_matrix:
            similarity_matrix = self.similarity_matrix(fps)

        self.clusters = pd.Series(self._cluster(similarity_matrix),
                                  index=fps.index,
                                  name='clusters')

    def _cluster_cumsum(self, shuffled=True):

        nums = self.clusters.value_counts()
        if shuffled:
            nums = nums.ix[np.random.permutation(nums.index)].cumsum()
        return nums

    def split(self, ratio):

        """ Return splits of the data with thresholded similarity according to a
        specified ratio.

        Args:
            ratio (tuple[ints]):
                the ratio to use.
        Returns:
            generator[pd.Series]:
                Generator of boolean split masks for the reqested splits.

        Example:
            st = SimThresholdSplit(ms, fper='morgan', similarity_metric='jaccard')
            train, valid, test = st.split(ratio=(70, 15, 15))
        """

        ratio = self.split_sizes(ratio)
        nums = self._cluster_cumsum()
        res = pd.Series(np.nan, index=nums.index, name='split')

        for i, _ in enumerate(ratio):
            lower = 0 if i == 0 else sum(ratio[:i])
            upper = len(ratio) if i == len(ratio) else sum(ratio[:i + 1])
            res[nums[(nums > lower) & (nums <= upper)].index] = i

        res = res.sort_index()
        res = self.clusters.to_frame().join(res, on='clusters')['split']
        return (res == i for i, _ in enumerate(ratio))

    def k_fold(self, n_folds):

        """ Returns k-fold cross-validated folds with thresholded similarity.

        Args:
            n_folds (int):
                The number of folds to provide.

        Returns:
            generator[(pd.Series, pd.Series)]:
                The splits in series.
        """

        folds = self.split((1,) * n_folds)
        return ((~fold, fold) for fold in folds)


    def split_sizes(self, ratio):
        """ Calculate the sizes of the splits """

        tot = sum(ratio)
        return [self.n_instances * rat / tot for rat in ratio]


    def similarity_matrix(self, fps):
        """ Calculate the similarity matrix for fingerprints. """

        if self.memory_optimized:
            return self._sim_low_mem(fps)
        else:
            return self._sim(fps)


    def _sim(self, fps):
        """ Fast but memory intensive implementation of similarity matrix
        calculation. """

        D = squareform(pdist(fps, self.similarity_metric))
        D = 1 - D # similarity is 1 - distance
        return triu(D >= self.threshold, k=1).todok()

    def _sim_low_mem(self, fps):
        """ Slow but memory efficient implementation of similarity matrix
        calculation """

        S = dok_matrix(len(fps), len(fps))
        for i, fp in enumerate(fps):
            D = cdist(fp[np.newaxis, :], fps[i + 1:], self.similarity_metric)
            D = 1 - D
            S[i, i + 1:] = dok_matrix(D >= threshold)
        return S

    def _cluster(self, S):
        """ Assign instances to clusters. """

        pairs = sorted(S.keys(), key=lambda x: x[0]) # sort pairs by first index
        clustered = np.arange(self.n_instances)

        for i, j in pairs:
            clustered[j] = clustered[i]

        return clustered
