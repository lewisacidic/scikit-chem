import pandas as pd
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
import gzip
try:
    import cPickle as pickle
except:
    import pickle
import skchem
import os
from skchem.target_prediction import TargetPredictionAlgorithm

class PIDGIN(object):
    
    def __init__(self):
        with gzip.open(os.path.join(os.path.dirname(__file__),'../data/PIDGIN_models.pkl.gz'), 'rb') as f:
            self.models = pickle.load(f)
        self.fingerprint = skchem.skchemize(GetMorganFingerprintAsBitVect, radius=2, nBits=2048)
        self.targets = self.models.keys()

    def __call__(self, m):
        return self.predict_proba(m)
        
    def predict(self, m):

        """ Predict binary binding profile for a molecule against 1080 protein targets """

        fp = self.fingerprint(m)
        res = pd.Series(index=self.targets)
        for target, model in self.models.iteritems():
            res[target] = model.predict(fp)[0]
        return res

    def map_predict(self, m):
        fp = self.fingerprint(m)
        return pd.Series(map(lambda k: self.models[k].predict(fp), self.targets), index=self.targets)

    def predict_proba(self, m):

        """ Predict probability of molecule m binding to 1080 protein targets """

        fp = self.fingerprint(m)
        res = pd.Series(index=self.targets)
        for target, model in self.models.iteritems():
            res[target] = model.predict_proba(fp)[:, 1][0]
        return res

    def map_predict_proba(self, m):

        """ Predict the log probability of molecule m binding to the 1080 proteins """

        fp = self.fingerprint(m)
        return pd.Series(map(lambda k: self.models[k].predict_proba(fp)[:, 1][0], self.targets), index=self.targets)

    def predict_log_proba(self, m):

        """ Predict the log probability of molecule m binding to the 1080 proteins """

        fp = self.fingerprint(m)
        res = pd.Series(index=self.targets)
        for target, model in self.models.iteritems():
            res[target] = model.predict_log_proba(fp)[:, 1][0]
        return res

    def map_predict_log_proba(self, m):

        """ Predict the log probabiltiy of molecule m binding to the 1080 proteins using map, for simple parallelism """

        fp = self.fingerprint(m)
        return pd.Series(map(lambda k: self.models[k].predict_log_proba(fp[:, 1][0]), self.targets), index=self.targets)

    def df_predict(self, df):

        """more efficient way to call the predict on large scikit-chem style dataframes"""

        fps = df.structure.apply(self.fingerprint)
        res = pd.DataFrame(index=fps.index, columns=self.targets)
        for target, model in self.models.iteritems():
            res[target] = model.predict(fps)
        return res

    def df_map_predict(self, df):

        """ More efficient way to call the predict on large scikit-chem style dataframes, with a map implementation for easy parallelism"""

        fps = df.structure.apply(self.fingerprint)

        return pd.DataFrame(map(lambda k: self.models[k].predict(fps), self.targets), columns=fps.index, index=self.targets).T


    def df_predict_proba(self, df):
        
        """ More efficient way to call the predict_proba on large scikit-chem style dataframes"""


        fps = df.structure.apply(self.fingerprint)
        res = pd.DataFrame(index=fps.index, columns=self.targets)

        #parallelize here
        for target, model in self.models.iteritems():
            res[target] = model.predict_proba(fps)[:, 1]
        return res

    def df_map_predict_proba(self, df):

        ''' map based way to call the predict_proba on large scikit-chem style dataframes'''

        fps = df.structure.apply(self.fingerprint)

        #parallize here trivially
        return pd.DataFrame(map(lambda k: self.models[k].predict_proba(fps)[:, 1], self.targets), columns=fps.index, index=self.targets).T

    def df_predict_log_proba(self, df):

        """ More efficient way to call the predict_proba on large scikit-chem style dataframes"""

        fps = df.structure.apply(self.fingerprint)
        res = pd.DataFrame(index=fps.index, columns=self.targets)

        for target, model in self.models.iteritems():
            res[target] = model.predict_log_proba(fps)[:, 1]
        return res

    def df_map_predict_log_proba(self, df):

        """ More efficient way to call the predict on large scikit-chem style dataframes, with a map implementation for easy parallelism"""

        fps = df.structure.apply(self.fingerprint)
        
        return pd.DataFrame(map(lambda k: self.models[k].predict_log_proba(fps)[:, 1], self.targets), columns=fps.index, index=self.targets).T