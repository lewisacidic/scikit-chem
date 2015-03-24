#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.target_prediction.target_prediction

Target prediction base class
"""

import pandas as pd
from skchem import Mol

class AbstractTargetPredictionAlgorithm(object):

    """abstract target prediction class, inherit from this for a target prediction algorithm"""

    def predict(self, obj):

        """ Predict the targets for a compound or compound DataFrame """

        if isinstance(obj, Mol):
            return self._m_predict(obj)

        if isinstance(obj, pd.Series) or isinstance(obj, pd.DataFrame):
            return self._df_predict(obj)

    def predict_proba(self, obj):

        """ Predict the probability of hitting targets for a compound or compound DataFrame """

        if isinstance(obj, Mol):
            return self._m_predict_proba(obj)

        if isinstance(obj, pd.Series) or isinstance(obj, pd.DataFrame):
            return self._df_predict_proba(obj)

    def predict_log_proba(self, obj):

        """ Predict the log probability of hitting target for a compound or compound DataFrame """

        if isinstance(obj, Mol):
            return self._m_predict_log_proba(obj)

        if isinstance(obj, pd.Series) or isinstance(obj, pd.DataFrame):
            return self._df_predict_log_proba(obj)

    def _m_predict(self, m):

        """ Overload """

        raise NotImplementedError

    def _m_predict_proba(self, m):

        """ Overload """

        raise NotImplementedError

    def _m_predict_log_proba(self, m):

        """ Overload """

        raise NotImplementedError

    def _df_predict(self, df):

        """ Overload for prediction """

        raise NotImplementedError

    def _df_predict_proba(self, df):

        """ Overload for prediction """

        raise NotImplementedError

    def _df_predict_log_proba(self, df):

        """ Overload for prediction """

        raise NotImplementedError

    def __call__(self, obj):

        """ Predict the targets for a compound or compound DataFrame """

        self.predict(obj)
