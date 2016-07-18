#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.atom

Module specifying atom based descriptor generators.
"""

import logging
import subprocess
import tempfile

import pandas as pd
import numpy as np
import time

from ..io import write_sdf
from .. import core
from ..utils import NamedProgressBar, line_count

LOGGER = logging.getLogger(__file__)

# TODO: fix averagemicrospeciescharge
# TODO: fix logd logp logs
# TODO: oen (orbital electronegativity) - sigma + pi
# TODO: water accessible surface area

class ChemAxonFeatureCalculator(object):
    _optimal_feats = ['acceptorcount', 'accsitecount', 'aliphaticatomcount', 'aliphaticbondcount', 'aliphaticringcount',
                      'aromaticatomcount', 'aromaticbondcount', 'aromaticringcount', 'asymmetricatomcount',
                      'averagemolecularpolarizability', 'axxpol', 'ayypol', 'azzpol', 'balabanindex', 'bondcount',
                      'carboaliphaticringcount', 'carboaromaticringcount', 'carboringcount', 'chainatomcount',
                      'chainbondcount', 'charge', 'chiralcentercount', 'connectedgraph', 'cyclomaticnumber', 'dipole',
                      'donorcount', 'donorsitecount', 'doublebondstereoisomercount', 'dreidingenergy', 'formalcharge',
                      'fragmentcount', 'fsp3', 'fusedaliphaticringcount', 'fusedaromaticringcount', 'fusedringcount',
                      'hararyindex', 'heteroaliphaticringcount', 'heteroaromaticringcount', 'heteroringcount', 'hlb',
                      'hmopienergy', 'hyperwienerindex', 'largestringsize', 'largestringsystemsize',
                      'markushenumerationcount', 'maximalprojectionarea', 'maximalprojectionradius',
                      'maximalprojectionsize', 'minimalprojectionarea', 'minimalprojectionradius',
                      'minimalprojectionsize', 'mmff94energy', 'molpol', 'pienergy', 'plattindex', 'psa', 'randicindex',
                      'refractivity', 'resonantcount', 'ringatomcount', 'ringbondcount', 'ringcount', 'ringsystemcount',
                      'rotatablebondcount', 'smallestatomringsize', 'smallestringsize', 'smallestringsystemsize',
                      'stereodoublebondcount', 'stereoisomercount', 'szegedindex', 'tetrahedralstereoisomercount',
                      'vdwsa', 'volume', 'wateraccessiblesurfacearea', 'wienerindex', 'wienerpolarity']

    _all_feats = ['acc', 'acceptor', 'acceptorcount', 'acceptormultiplicity', 'acceptorsitecount', 'acceptortable',
                  'accsitecount', 'aliphaticatom', 'aliphaticatomcount', 'aliphaticbondcount', 'aliphaticringcount',
                  'aliphaticringcountofsize', 'aromaticatom', 'aromaticatomcount', 'aromaticbondcount',
                  'aromaticelectrophilicityorder', 'aromaticnucleophilicityorder', 'aromaticringcount',
                  'aromaticringcountofsize', 'asa', 'asymmetricatom', 'asymmetricatomcount', 'asymmetricatoms',
                  'atomicpolarizability', 'atompol', 'averagemicrospeciescharge', 'averagemolecularpolarizability',
                  'averagepol', 'avgpol', 'axxpol', 'ayypol', 'azzpol', 'balabanindex', 'bondcount',
                  'canonicalresonant', 'canonicaltautomer', 'carboaliphaticringcount', 'carboaromaticringcount',
                  'carboringcount', 'chainatom', 'chainatomcount', 'chainbondcount', 'charge', 'chargedensity',
                  'chargedistribution', 'chiralcenter', 'chiralcentercount', 'chiralcenters', 'connectedgraph',
                  'cyclomaticnumber', 'dipole', 'distancedegree', 'don', 'donor', 'donorcount', 'donormultiplicity',
                  'donorsitecount', 'donortable', 'donsitecount', 'doublebondstereoisomercount', 'dreidingenergy',
                  'eccentricity', 'electrondensity', 'electrophilicityorder', 'electrophiliclocalizationenergy',
                  'enumerationcount', 'enumerations', 'formalcharge', 'fragmentcount', 'fsp3',
                  'fusedaliphaticringcount', 'fusedaromaticringcount', 'fusedringcount', 'generictautomer',
                  'hararyindex', 'hasvalidconformer', 'hbda', 'hbonddonoracceptor', 'heteroaliphaticringcount',
                  'heteroaromaticringcount', 'heteroringcount', 'hindrance', 'hlb', 'hmochargedensity',
                  'hmoelectrondensity', 'hmoelectrophilicityorder', 'hmoelectrophiliclocalizationenergy', 'hmohuckel',
                  'hmohuckeleigenvalue', 'hmohuckeleigenvector', 'hmohuckelorbitals', 'hmohuckeltable',
                  'hmolocalizationenergy', 'hmonucleophilicityorder', 'hmonucleophiliclocalizationenergy',
                  'hmopienergy', 'huckel', 'huckeleigenvalue', 'huckeleigenvector', 'huckelorbitals', 'huckeltable',
                  'hyperwienerindex', 'ioncharge', 'isoelectricpoint', 'largestatomringsize', 'largestringsize',
                  'largestringsystemsize', 'localizationenergy', 'logd', 'logdcalculator', 'logp', 'logpcalculator',
                  'logs', 'majormicrospecies', 'majorms', 'majorms2', 'majortautomer', 'markushenumerationcount',
                  'markushenumerations', 'maximalprojectionarea', 'maximalprojectionradius', 'maximalprojectionsize',
                  'minimalprojectionarea', 'minimalprojectionradius', 'minimalprojectionsize', 'mmff94energy',
                  'molecularpolarizability', 'molecularsurfacearea', 'molpol', 'moststabletautomer', 'msa', 'msacc',
                  'msdon', 'name', 'nucleophilicityorder', 'nucleophiliclocalizationenergy', 'oen',
                  'orbitalelectronegativity', 'pi', 'pichargedensity', 'pienergy', 'pka', 'pkacalculator', 'pkat',
                  'plattindex', 'pol', 'polarizability', 'polarsurfacearea', 'psa', 'randicindex',
                  'randommarkushenumerations', 'refractivity', 'resonantcount', 'resonants', 'ringatom',
                  'ringatomcount', 'ringbondcount', 'ringcount', 'ringcountofatom', 'ringcountofsize',
                  'ringsystemcount', 'ringsystemcountofsize', 'rotatablebondcount', 'smallestatomringsize',
                  'smallestringsize', 'smallestringsystemsize', 'stereodoublebondcount', 'stereoisomercount',
                  'stericeffectindex', 'sterichindrance', 'szegedindex', 'tautomercount', 'tautomers',
                  'tetrahedralstereoisomercount', 'tholepolarizability', 'topanal', 'topologyanalysistable',
                  'totalchargedensity', 'tpol', 'tpolarizability', 'vdwsa', 'volume', 'wateraccessiblesurfacearea',
                  'wienerindex', 'wienerpolarity']

    def __init__(self, feat_set='optimal'):
        if feat_set == 'all':
            self.index = self._all_feats
        elif feat_set == 'optimal':
            self.index = self._optimal_feats
        elif feat_set in self._all_feats:
            self.index = [feat_set]
        elif isinstance(feat_set, (list, tuple)):
            valid = np.array([feat in self._all_feats for feat in feat_set])
            if all(valid):
                self.index = feat_set
            else:
                self.index = feat_set
                raise NotImplementedError('Descriptor \'{}\' not available.'.format(np.array(feat_set)[~valid]))
        else:
            raise NotImplementedError('Feature set {} not available.'.format(feat_set))

    def transform(self, obj):
        if isinstance(obj, core.Mol):
            return self._transform_series(pd.Series(obj)).iloc[0]
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_series(self, series):

        with tempfile.NamedTemporaryFile(suffix='.sdf') as in_file, tempfile.NamedTemporaryFile() as out_file:
            # write mols to file
            write_sdf(series, in_file.name)
            args = ['cxcalc', in_file.name, '-o', out_file.name] + self.index

            LOGGER.debug('Running: ' + ' '.join(args))

            # call command line
            bar = NamedProgressBar(name=self.__class__.__name__, max_value=len(series))
            p = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # monitor
            while p.poll() is None:
                bar.update(line_count(out_file.name))
                time.sleep(1)
            bar.update(len(series))
            p.wait()

            try:
                finished = pd.read_table(out_file.name).set_index('id')
            except Exception:
                finished = None
        finished.index = series.index
        return finished

class ChemAxonAtomFeatureCalculator(object):

    _all_feats = ['acceptormultiplicity', 'aliphaticatom', 'aromaticatom', 'aromaticelectrophilicityorder',
                  'asymmetricatom', 'atomicpolarizability', 'chainatom', 'chargedensity', 'chiralcenter',
                  'distancedegree', 'donormultiplicity', 'eccentricity', 'electrondensity',
                  'electrophiliclocalizationenergy', 'hindrance', 'hmochargedensity', 'hmoelectrondensity', 'hmoelectrophilicityorder',
                  'hmoelectrophiliclocalizationenergy', 'hmonucleophilicityorder', 'hmonucleophiliclocalizationenergy', 'ioncharge', 'largestatomringsize',
                  'nucleophilicityorder', 'nucleophiliclocalizationenergy', 'oen', 'pichargedensity', 'ringatom', 'ringcountofatom',
                  'stericeffectindex', 'totalchargedensity']

    _h_inc_feats = ['acc', 'atomicpolarizability', 'charge', 'distancedegree', 'don',
       'eccentricity', 'hindrance', 'largestatomringsize', 'oen',
       'ringcountofatom', 'smallestatomringsize', 'stericeffectindex']

    def __init__(self, feat_set='all', include_hs=False, max_atoms=75):

        """
        Args:
            feat_set (str or list<str>):
                The feature sets to calculate.
                - a single identifier as a `str`
                - a list of identifiers
                - 'h_inc' or those that also calculate for Hs.
                - 'all' for all

            max_atoms:
                - The maximum number of atoms available.
        """
        self.max_atoms = max_atoms

        if feat_set in self._all_feats:
            self.index = [feat_set]
        elif feat_set == 'h_inclusive':
            self.index = self._h_inc_feats
        elif feat_set == 'all':
            self.index = self._all_feats
        elif isinstance(feat_set, (list, tuple)):
            valid = np.array([feat in self._all_feats for feat in feat_set])
            if all(valid):
                self.index = feat_set
            else:
                raise NotImplementedError('Descriptor \'{}\' not available.'.format(np.array(feat_set)[~valid]))
            self.feature_names = feat_set
        else:
            raise NotImplementedError('{} feature set is not available'.format(feat_set))

    @property
    def feature_names(self):
        return self.index

    def transform(self, obj):
        if isinstance(obj, core.Atom):
            return self._transform_atom(obj)
        elif isinstance(obj, core.Mol):
            return self._transform_mol(obj)
        elif isinstance(obj, pd.Series):
            return self._transform_series(obj)
        elif isinstance(obj, pd.DataFrame):
            return self._transform_series(obj.structure)
        elif isinstance(obj, (tuple, list)):
            return self._transform_series(obj)
        else:
            raise NotImplementedError

    def _transform_atom(self, atom):
        raise NotImplementedError('Cannot calculate atom wise with Chemaxon')

    def _transform_mol(self, mol):
        # make into series then use self._transform_mol
        ser = pd.Series([mol], name=mol.name)
        res = self._transform_series(ser)
        return res.iloc[0]

    def _transform_series(self, series):

        with tempfile.NamedTemporaryFile(suffix='.sdf') as in_file, tempfile.NamedTemporaryFile() as out_file:
            # write mols to file
            write_sdf(series, in_file.name)
            args = ['cxcalc', in_file.name, '-o', out_file.name] + self.index

            LOGGER.debug('Running: ' + ' '.join(args))

            # call command line
            p = subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            bar = NamedProgressBar(name=self.__class__.__name__, max_value=len(series))

            while p.poll() is None:
                bar.update(line_count(out_file.name))
                time.sleep(1)
            bar.update(len(series))
            p.wait()
            finished = pd.read_table(out_file.name).set_index('id')

        def to_padded(s):
            res = np.repeat(np.nan, self.max_atoms)

            def parse_string(s):
                if s == '':
                    return np.nan
                elif s == 'false':
                    return 0
                elif s == 'true':
                    return 1
                else:
                    return float(s)

            ans = np.array([parse_string(i) for i in s.split(';')])
            res[:len(ans)] = ans
            return res
        res = np.array([[to_padded(i) for k, i in val.items()] for idx, val in finished.T.items()])
        res = pd.Panel(res, items=series.index, major_axis=pd.Index(finished.columns, name='cx_atom_desc'),
                        minor_axis=pd.Index(range(self.max_atoms), name='atom_idx'))
        return res.swapaxes(1, 2) # to be consistent with AtomFeatureCalculator
