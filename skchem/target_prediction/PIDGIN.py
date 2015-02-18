import pandas as pd
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
import gzip
import cPickle
import skchem
import os

class PIDGIN(object):
    
    def __init__(self):
        with gzip.open(os.path.join(os.path.dirname(__file__),'../data/PIDGIN_models.pkl.gz'), 'rb') as f:
            self.models = cPickle.load(f)
        self.fingerprint = skchem.skchemize(GetMorganFingerprintAsBitVect, 2, nBits=2048)

    def __call__(self, m):
        return self.predict_proba(m)
        
    def predict(self, m):
        fp = self.fingerprint(m)
        res = pd.Series()
        for target, model in self.models.iteritems():
            res[target] = model.predict(fp)[0]
        return res

    def predict_proba(self, m):

        '''predict probability of molecule m binding to 1080 proteins'''

        fp = self.fingerprint(m)
        res = pd.Series()
        for target, model in self.models.iteritems():
            res[target] = model.predict_proba(fp)[:,1][0]
        return res

    def df_predict(self, df):

        '''more efficient way to call the predict on large scikit-chem style dataframes'''

        fps = df.structure.apply(self.fingerprint)
        res = pd.DataFrame(index=fps.index, columns=self.models.keys())
        for target, model in self.models.iteritems():
            res[target] = model.predict(fps)
        return res

    def df_predict_proba(self, df):
        
        '''more efficient way to call the predict_proba on large scikit-chem style dataframes'''

        fps = df.structure.apply(self.fingerprint)
        res = pd.DataFrame(index=fps.index, columns=self.models.keys())

        #parallelize here
        for target, model in self.models.iteritems():
            res[target] = model.predict_proba(fps)[:, 1]
        return res

    def df_map_predict_proba(self, df):

        ''' map based way to call the predict_proba on large scikit-chem style dataframes'''

        fps = df.structure.apply(self.fingerprint)
        ts = self.models.keys()

        #parallize here trivially
        return pd.DataFrame(map(lambda k: self.models[k].predict_proba(fps)[:, 1], ts), columns=fps.index, index=ts).T


if __name__ == "__main__":

    from rdkit import Chem
    from time import time

    p = PIDGIN()
    m = Chem.MolFromSmiles('c1ccccc1')

    t = time()
    print p.predict_proba(m)
    print 'Prediction for one molecule took', time() - t, 'seconds'

    df = skchem.read_sdf('/Users/RichLewis/Dropbox/PhD/Data/sa_named.sdf')
    t = time()
    print p.df_predict_proba(df)
    print 'Prediction for', df.shape[0], 'molecules took', time() - t, 'seconds'
    
    t = time()
    print p.df_map_predict_proba(df)
    print 'Prediction for', df.shape[0], 'molecules using map took', time() - t, 'seconds'

    t = time()
    print df.structure.apply(p.predict_proba)
    print 'Prediction for', df.shape[0], 'molecules using slow method took', time() - t, 'seconds'