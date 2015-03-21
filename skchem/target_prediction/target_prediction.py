#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import pandas as _pd
from skchem import Mol as _Mol

class TargetPredictionAlgorithm(object):

    """abstract target prediction class, inherit from this for a target prediction algorithm"""

    def predict(self, obj):
        if type(obj) is _Mol:
            return self._m_predict(obj)

        if type(obj) in [_pd.Series, _pd.DataFrame]:
            return self._df_predict(obj)

    def predict_proba(self, obj):
        if type(obj) is _Mol:
            return self._m_predict_proba(obj)

        if type(obj) in [_pd.Series, _pd.DataFrame]:
            return self._df_predict_proba(obj)

    def predict_log_proba(self, obj):
        if type(obj) is _Mol:
            return self._m_predict_log_proba(obj)

        if type(obj) in [_pd.Series, _pd.DataFrame]:
            return self._df_predict_log_proba(obj)

    def _m_predict(self, m):
        return _pd.Series([0])

    def _m_predict_proba(self, m):
        return _pd.Series([0])

    def _m_predict_log_proba(self, m):
        return _pd.Series([0])

    def _df_predict(self, df):
        return _pd.Series([0])

    def _df_predict_proba(self, df):
        return _pd.Series([0])

    def _df_predict_proba(self, df):
        return _pd.Series([0])

    def __call__(self, obj):
        self.predict(obj)