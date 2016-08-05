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
import re
from abc import ABCMeta

import pandas as pd
import numpy as np

from ..utils import line_count, nanarray
from ..base import CLIWrapper, Transformer, AtomTransformer, BatchTransformer, Featurizer

LOGGER = logging.getLogger(__file__)

# TODO: fix averagemicrospeciescharge
# TODO: fix logd logp logs
# TODO: oen (orbital electronegativity) - sigma + pi
# TODO: water accessible surface area

#  TODO: these don't produce csv
# ['doublebondstereoisomers', 'conformers', 'stereoisomers', 'moleculardynamics', 'stereoanalysis',
# 'lowestenergyconformer', 'msdistr2', 'conformations', 'dominanttautomerdistribution', 'hnmr', 'moldyn', 'cnmr',
# 'frameworks', 'microspeciesdistribution', 'nmr', 'leconformer', 'msdistr', 'tetrahedralstereoisomers']
#

CHEMAXON_HINT = """ Install ChemAxon from https://www.chemaxon.com.  It requires a license, which can be freely obtained
for academics. """


class ChemAxonBaseFeaturizer(CLIWrapper, Featurizer):

    __metaclass__ = ABCMeta

    install_hint = CHEMAXON_HINT

    _feat_columns = {'averagepol': ['a_avg'], 'name': ['preferred_iupac_name'],
                     'aromaticbondcount': ['aromatic_bond_count'],
                     'maximalprojectionradius': ['maximal_projection_radius'],
                     'tpolarizability': ['a_avg', 'a_xx', 'a_yy', 'a_zz'], 'distance': ['distance'],
                     'acceptor': ['acceptor_count', 'acceptor_site_count'], 'fusedringcount': ['fused_ring_count'],
                     'charge': ['total_charge'], 'donor': ['donor_count', 'donor_site_count'], 'ringcount': ['ring_count'],
                     'chainbond': ['chain_bond'], 'mmff94energy': ['mmff94_energy'],
                     'huckel': ['aromatic_e+/nu-_order', 'localization_energy_l_+/l-', 'pi_energy',
                                'electron_density', 'charge_density'], 'chainatom': ['chain_atom'],
                     'shortestpath': ['shortest_path'], 'resonantcount': ['resonant_count'],
                     'tpol': ['a_avg', 'a_xx', 'a_yy', 'a_zz'], 'moststabletautomer': ['most_stable_tautomer'],
                     'generictautomer': ['generictautomer'], 'hmoelectrophilicityorder': ['hmo_aromatic_e+_order'],
                     'ringsystemcountofsize': ['ring_system_count_of_size'],
                     'largestatomringsize': ['largest_ring_size_of_atom'],
                     'tetrahedralstereoisomercount': ['tetrahedral_stereoisomer_count'], 'enumerations': ['enumerations'],
                     'ringatom': ['ring_atom'], 'connected': ['connected'],
                     'hmolocalizationenergy': ['hmo_localization_energy_l+/l-'],
                     'averagemolecularpolarizability': ['a_avg'], 'donorsitecount': ['donor_site_count'],
                     'donorcount': ['donor_count'], 'asymmetricatom': ['asymmetric_atom'], 'pienergy': ['pi_energy'],
                     'bondcount': ['bond_count'], 'chiralcenters': ['chiral_centers'],
                     'hmohuckel': ['hmo_aromatic_e+/nu-_order', 'hmo_localization_energy_l+/l-', 'hmo_pi_energy',
                                   'hmo_electron_density', 'hmo_charge_density'], 'huckeleigenvector': ['eigenvector'],
                     'ringcountofsize': ['ring_count_of_size'],
                     'heteroaliphaticringcount': ['heteroaliphatic_ring_count'], 'markushenumerations': ['enumerations'],
                     'minimalprojectionradius': ['minimal_projection_radius'], 'dipole': ['dipole'],
                     'balabanindex': ['balaban_index'], 'aromaticnucleophilicityorder': ['aromatic_nu-_order'],
                     'tautomercount': ['tautomer_count'], 'cyclomaticnumber': ['cyclomatic_number'],
                     'psa': ['polar_surface_area'], 'isoelectricpoint': ['pi'], 'hmopienergy': ['hmo_pi_energy'],
                     'ayypol': ['a_yy'], 'fragmentcount': ['fragment_count'], 'acceptormultiplicity': ['acceptor_multiplicity'],
                     'topologyanalysistable': ['atom_count', 'aliphatic_atom_count', 'aromatic_atom_count',
                                               'bond_count', 'aliphatic_bond_count', 'aromatic_bond_count',
                                               'rotatable_bond_count', 'ring_count', 'aliphatic_ring_count',
                                               'aromatic_ring_count', 'hetero_ring_count', 'heteroaliphatic_ring_count',
                                               'heteroaromatic_ring_count', 'ring_atom_count', 'ring_bond_count',
                                               'chain_atom_count', 'chain_bond_count', 'smallest_ring_size',
                                               'largest_ring_size'], 'ioncharge': ['charge'],
                     'asymmetricatoms': ['asymmetric_atoms'],
                     'wateraccessiblesurfacearea': ['asa', 'asa+', 'asa-', 'asa_h', 'asa_p'], 'avgpol': ['a_avg'],
                     'carboaliphaticringcount': ['carboaliphatic_ring_count'],
                     'aliphaticringcount': ['aliphatic_ring_count'], 'donormultiplicity': ['donor_multiplicity'],
                     'minimalprojectionarea': ['minimal_projection_area'],
                     'nucleophiliclocalizationenergy': ['localization_energy_l-'], 'dihedral': ['dihedral'],
                     'heteroringcount': ['hetero_ring_count'], 'azzpol': ['a_zz'],
                     'molecularsurfacearea': ['van_der_waals_surface_area_3d'],
                     'hmonucleophiliclocalizationenergy': ['hmo_localization_energy_l-'],
                     'chargedistribution': ['charge_distribution'], 'pol': ['a_mol', 'a_atom'],
                     'hmoelectrondensity': ['hmo_electron_density'], 'carboaromaticringcount': ['carboaromatic_ring_count'],
                     'acceptorsitecount': ['acceptor_site_count'], 'markushenumerationcount': ['markush_library_size'],
                     'localizationenergy': ['localization_energy_l+/l-'], 'hararyindex': ['harary_index'],
                     'asa': ['asa', 'asa+', 'asa-', 'asa_h', 'asa_p'], 'acc': ['acc'], 'majortautomer': ['major_tautomer'],
                     'majormicrospecies': ['major-ms'], 'aliphaticatomcount': ['aliphatic_atom_count'],
                     'angle': ['angle'], 'huckeleigenvalue': ['eigenvalue'], 'axxpol': ['a_xx'],
                     'chiralcenter': ['chiral_center'], 'aliphaticbondcount': ['aliphatic_bond_count'],
                     'smallestatomringsize': ['smallest_ring_size_of_atom'], 'dreidingenergy': ['dreiding_energy'],
                     'maximalprojectionsize': ['length_perpendicular_to_the_max_area'],
                     'largestringsystemsize': ['largest_ring_system_size'], 'accsitecount': ['acceptor_site_count'],
                     'refractivity': ['refractivity'], 'bondtype': ['bond_type'], 'chargedensity': ['charge_density'],
                     'resonants': ['resonants'], 'aromaticatomcount': ['aromatic_atom_count'],
                     'distancedegree': ['distance_degree'], 'hasvalidconformer': ['has_valid_conformer'],
                     'electrondensity': ['electron_density'], 'asymmetricatomcount': ['asymmetric_atom_count'],
                     'fsp3': ['fsp3'], 'don': ['don'], 'fusedaliphaticringcount': ['fused_aliphatic_ring_count'],
                     'pkat': ['pkat'], 'fusedaromaticringcount': ['fused_aromatic_ring_count'],
                     'majorms2': ['majorms2'], 'maximalprojectionarea': ['maximal_projection_area'],
                     'hbonddonoracceptor': ['acceptor_count', 'donor_count', 'acceptor_site_count', 'donor_site_count'],
                     'acceptorcount': ['acceptor_count'], 'molecularpolarizability': ['a_mol'],
                     'huckeltable': ['aromatic_e+/nu-_order', 'localization_energy_l+/l-', 'pi_energy',
                                     'electron_density', 'charge_density'],
                     'rotatablebondcount': ['rotatable_bond_count'],
                     'minimalprojectionsize': ['length_perpendicular_to_the_min_area'],
                     'polarizability': ['a_mol', 'a_atom'], 'acceptortable': ['acceptor_count', 'acceptor_site_count'],
                     'aliphaticringcountofsize': ['aliphatic_ring_count_of_size'], 'hlb': ['hlb'],
                     'eccentricity': ['eccentricity'], 'hmochargedensity': ['hmo_charge_density'],
                     'hmohuckeleigenvalue': ['hmo_eigenvalue'], 'totalchargedensity': ['total_charge_density'],
                     'hmonucleophilicityorder': ['hmo_aromatic_nu-_order'],
                     'aromaticringcountofsize': ['aromatic_ring_count_of_size'],
                     'electrophilicityorder': ['aromatic_e+_order'], 'connectedgraph': ['connected_graph'],
                     'plattindex': ['platt_index'], 'logp': ['logp'],
                     'topanal': ['atom_count', 'aliphatic_atom_count', 'aromatic_atom_count', 'bond_count',
                                 'aliphatic_bond_count', 'aromatic_bond_count', 'rotatable_bond_count', 'ring_count',
                                 'aliphatic_ring_count', 'aromatic_ring_count', 'hetero_ring_count',
                                 'heteroaliphatic_ring_count', 'heteroaromatic_ring_count', 'ring_atom_count',
                                 'ring_bond_count', 'chain_atom_count', 'chain_bond_count', 'smallest_ring_size',
                                 'largest_ring_size'],
                     'logdcalculator': ['ph=0', 'ph=1', 'ph=2', 'ph=3', 'ph=4', 'ph=5', 'ph=6', 'ph=7', 'ph=8', 'ph=9',
                                        'ph=10', 'ph=11', 'ph=12', 'ph=13', 'ph=14', 'unnamed:_16'],
                     'logs': ['ph=0.0', 'ph=1.0', 'ph=2.0', 'ph=3.0', 'ph=4.0', 'ph=5.0', 'ph=6.0', 'ph=7.0', 'ph=8.0',
                              'ph=9.0', 'ph=10.0', 'ph=11.0', 'ph=12.0', 'ph=13.0', 'ph=14.0', 'unnamed:_16'],
                     'atompol': ['a_atom'], 'canonicalresonant': ['structure'], 'ringbond': ['ring_bond'],
                     'ringatomcount': ['ring_atom_count'], 'donortable': ['donor_count', 'donor_site_count'],
                     'randicindex': ['randic_index'], 'rotatablebond': ['rotatable_bond'],
                     'hyperwienerindex': ['hyper_wiener_index'], 'hmohuckeleigenvector': ['hmo_eigenvector'],
                     'carboringcount': ['carbo_ring_count'], 'logpcalculator': ['logp', 'unnamed:_2'],
                     'ringsystemcount': ['ring_system_count'], 'largestringsize': ['largest_ring_size'],
                     'stereodoublebondcount': ['stereo_double_bond_count'], 'pi': ['pi'],
                     'stericeffectindex': ['steric_effect_index'], 'volume': ['van_der_waals_volume'],
                     'averagemicrospeciescharge': ['charge'], 'pka': ['apka1', 'apka2', 'bpka1', 'bpka2', 'atoms'],
                     'hmohuckeltable': ['hmo_aromatic_e+/nu-_order', 'hmo_localization_energy_l+/l-', 'hmo_pi_energy',
                                        'hmo_electron_density', 'hmo_charge_density'],
                     'ringcountofatom': ['ring_count_of_atom'],
                     'aromaticelectrophilicityorder': ['aromatic_e+_order'], 'hindrance': ['steric_hindrance'],
                     'chainatomcount': ['chain_atom_count'],
                     'pkacalculator': ['apka1', 'apka2', 'bpka1', 'bpka2', 'atoms'],
                     'heteroaromaticringcount': ['heteroaromatic_ring_count'], 'sterichindrance': ['steric_hindrance'],
                     'hbda': ['acceptor_count', 'donor_count', 'acceptor_site_count', 'donor_site_count'], 'molpol': ['a_mol'],
                     'atomicpolarizability': ['a_atom'],
                     'msdon': ['ph=0.00', 'ph=1.00', 'ph=2.00', 'ph=3.00', 'ph=4.00', 'ph=5.00', 'ph=6.00', 'ph=7.00',
                               'ph=8.00', 'ph=9.00', 'ph=10.00', 'ph=11.00', 'ph=12.00', 'ph=13.00', 'ph=14.00'],
                     'enumerationcount': ['markush_library_size'], 'vdwsa': ['van_der_waals_surface_area_3d'],
                     'orbitalelectronegativity': ['sigma_orbital_electronegativity', 'pi_orbital_electronegativity'],
                     'hmoelectrophiliclocalizationenergy': ['hmo_localization_energy_l+'],
                     'smallestringsize': ['smallest_ring_size'], 'szegedindex': ['szeged_index'],
                     'nucleophilicityorder': ['aromatic_nu-_order'], 'canonicaltautomer': ['canonical_tautomer'],
                     'stereoisomercount': ['stereoisomer_count'], 'msa': ['van_der_waals_surface_area_3d'],
                     'donsitecount': ['donor_site_count'], 'randommarkushenumerations': ['randommarkushenumerations'],
                     'wienerindex': ['wiener_index'], 'huckelorbitals': ['orbitals'],
                     'doublebondstereoisomercount': ['double_bond_stereoisomer_count'], 'tautomers': ['tautomers'],
                     'polarsurfacearea': ['polar_surface_area'], 'chiralcentercount': ['chiral_center_count'],
                     'electrophiliclocalizationenergy': ['localization_energy_l+'],
                     'aliphaticatom': ['aliphatic_atom'], 'ringbondcount': ['ring_bond_count'],
                     'wienerpolarity': ['wiener_polarity'],
                     'msacc': ['ph=0.00', 'ph=1.00', 'ph=2.00', 'ph=3.00', 'ph=4.00', 'ph=5.00', 'ph=6.00', 'ph=7.00',
                               'ph=8.00', 'ph=9.00', 'ph=10.00', 'ph=11.00', 'ph=12.00', 'ph=13.00', 'ph=14.00'],
                     'formalcharge': ['formal_charge'], 'smallestringsystemsize': ['smallest_ring_system_size'],
                     'majorms': ['major-ms'], 'tholepolarizability': ['a_avg', 'a_xx', 'a_yy', 'a_zz'],
                     'aromaticatom': ['aromatic_atom'],
                     'oen': ['sigma_orbital_electronegativity', 'pi_orbital_electronegativity'],
                     'chainbondcount': ['chain_bond_count'],
                     'logd': ['ph=0.00', 'ph=1.00', 'ph=2.00', 'ph=3.00', 'ph=4.00', 'ph=5.00', 'ph=6.00', 'ph=7.00',
                              'ph=8.00', 'ph=9.00', 'ph=10.00', 'ph=11.00', 'ph=12.00', 'ph=13.00', 'ph=14.00'],
                     'hmohuckelorbitals': ['hmo_orbitals'], 'aromaticringcount': ['aromatic_ring_count'],
                     'pichargedensity': ['pi_charge_density']}

    def __init__(self, features='optimal', **kwargs):
        super(ChemAxonBaseFeaturizer, self).__init__(**kwargs)
        self.features = features

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, features):
        if features in ('optimal', 'all'):
            self._features = self._optimal_feats
        elif isinstance(features, str):
            self.features = [features]
        elif isinstance(features, (list, tuple)):
            valid = np.array([feat in self._feat_columns.keys() for feat in features])
            if not all(valid):
                raise NotImplementedError('Descriptor \'{}\' not available.'.format(np.array(features)[~valid]))
            else:
                self._features = list(features)
        else:
            raise NotImplementedError('Feature set {} not available.'.format(features))

    def _feature_index(self):
        return pd.Index(sum((self._feat_columns[feat] for feat in self.features), []), name='features')

    def validate_install(self):
        try:
            return 0 == subprocess.call(['cxcalc'], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        except FileNotFoundError:
            return False

    def monitor_progress(self, filename):
        res = line_count(filename)
        return res - 1 if res else 0

    def _cli_args(self, infile, outfile):
        return ['cxcalc', infile, '-o', outfile] + self.features

    def _parse_outfile(self, outfile):
        res = pd.read_table(outfile, engine='python').drop('id', axis=1)
        return res

    def _parse_errors(self, errs):
        LOGGER.debug('stderr: %s', errs)
        return [] # instances are not skipped ever, so don't return anything


class ChemAxonFeaturizer(ChemAxonBaseFeaturizer, BatchTransformer, Transformer):

    _optimal_feats = ['acceptorcount', 'accsitecount', 'aliphaticatomcount', 'aliphaticbondcount', 'aliphaticringcount',
                      'aromaticatomcount', 'aromaticbondcount', 'aromaticringcount', 'asymmetricatomcount',
                      'averagemolecularpolarizability', 'axxpol', 'ayypol', 'azzpol', 'balabanindex', 'bondcount',
                      'carboaliphaticringcount', 'carboaromaticringcount', 'carboringcount', 'chainatomcount',
                      'chainbondcount', 'chiralcentercount', 'connectedgraph', 'cyclomaticnumber', 'dipole',
                      'donorcount', 'donorsitecount', 'doublebondstereoisomercount', 'dreidingenergy', 'formalcharge',
                      'fragmentcount', 'fsp3', 'fusedaliphaticringcount', 'fusedaromaticringcount', 'fusedringcount',
                      'hararyindex', 'heteroaliphaticringcount', 'heteroaromaticringcount', 'heteroringcount', 'hlb',
                      'hmopienergy', 'hyperwienerindex', 'largestringsize', 'largestringsystemsize',
                      'markushenumerationcount', 'maximalprojectionarea', 'maximalprojectionradius',
                      'maximalprojectionsize', 'minimalprojectionarea', 'minimalprojectionradius',
                      'minimalprojectionsize', 'mmff94energy', 'molpol', 'pienergy', 'plattindex', 'psa', 'randicindex',
                      'refractivity', 'resonantcount', 'ringatomcount', 'ringbondcount', 'ringcount', 'ringsystemcount',
                      'rotatablebondcount', 'smallestringsize', 'smallestringsystemsize',
                      'stereodoublebondcount', 'stereoisomercount', 'szegedindex', 'tetrahedralstereoisomercount',
                      'vdwsa', 'volume', 'wateraccessiblesurfacearea', 'wienerindex', 'wienerpolarity']

    @property
    def name(self):
        return 'cx_mol'

    @property
    def columns(self):
        return self._feature_index()

    def _parse_outfile(self, outfile):
        res = super(ChemAxonFeaturizer, self)._parse_outfile(outfile)
        return res.applymap(lambda s: np.nan if isinstance(s, str) and 'FAILED' in s else float(s))


class ChemAxonAtomFeaturizer(ChemAxonBaseFeaturizer, AtomTransformer, BatchTransformer):

    _optimal_feats = ['acceptormultiplicity', 'aliphaticatom', 'aromaticatom', 'aromaticelectrophilicityorder',
                      'asymmetricatom', 'atomicpolarizability', 'chainatom', 'chargedensity', 'chiralcenter',
                      'distancedegree', 'donormultiplicity', 'eccentricity', 'electrondensity',
                      'electrophiliclocalizationenergy', 'hindrance', 'hmochargedensity', 'hmoelectrondensity',
                      'hmoelectrophilicityorder',
                      'hmoelectrophiliclocalizationenergy', 'hmonucleophilicityorder',
                      'hmonucleophiliclocalizationenergy', 'ioncharge', 'largestatomringsize',
                      'nucleophilicityorder', 'nucleophiliclocalizationenergy', 'oen', 'pichargedensity', 'ringatom',
                      'ringcountofatom',
                      'stericeffectindex', 'totalchargedensity']

    _h_inc_feats = ['acc', 'atomicpolarizability', 'charge', 'distancedegree', 'don',
       'eccentricity', 'hindrance', 'largestatomringsize', 'oen',
       'ringcountofatom', 'smallestatomringsize', 'stericeffectindex']

    @property
    def name(self):
        return 'cx_atom'

    @property
    def minor_axis(self):
        return self._feature_index()

    def _transform_atom(self, atom):
        raise NotImplementedError('Cannot calculate per atom with ChemAxon')

    def _parse_outfile(self, outfile):
        res = super(ChemAxonAtomFeaturizer, self)._parse_outfile(outfile)

        def parse_string(s):
            if s == '':
                return np.nan
            elif s == 'false':
                return 0
            elif s == 'true':
                return 1
            else:
                try:
                    return float(s)
                except ValueError:
                    return np.nan

        def to_padded(s):
            res = np.repeat(np.nan, self.max_atoms)
            ans = np.array([parse_string(i) for i in str(s).split(';')])
            res[:len(ans)] = ans
            return res

        res = res.applymap(to_padded)
        return pd.Panel(res.values.tolist()).swapaxes(1, 2)


class ChemAxonNMRPredictor(ChemAxonBaseFeaturizer, BatchTransformer, AtomTransformer):

    _feat_columns = {'cnmr': ['cnmr'], 'hnmr': ['hnmr']}
    _optimal_feats = ['cnmr']

    def name(self):
        return 'cx_nmr'

    def _transform_atom(self, atom):
        raise NotImplementedError('ChemAxon cannot predict for atoms.')

    def monitor_progress(self, filename):
        return sum(1 for l in open(filename, 'rb') if l == b'##PEAKASSIGNMENTS=(XYMA)\r\n')

    @property
    def minor_axis(self):
        return pd.Index(self.features, name='shift')

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, val):
        if val == 'c':
            self._features = ['cnmr']
        elif val == 'h':
            self._features = ['hnmr']
        else:
            raise NotImplementedError('Feature {} not implemented'.format(val))

    def _parse_outfile(self, outfile):
        n_mols = self.monitor_progress(outfile)
        res = nanarray((n_mols, self.max_atoms, 1))
        regex = re.compile(b'\((-?\d+.\d+),\d+,[A-Z],<([0-9\,]+)>\)\r\n')

        mol_idx = 0

        with open(outfile, 'rb') as f:
            # loop through the file - inner loop will also advance the pointer
            for l in f:
                if l == b'##PEAKASSIGNMENTS=(XYMA)\r\n':
                    for row in f:
                        if row == b'##END=\r\n':
                            break
                        else:
                            LOGGER.debug('Row to parse: %s', row)
                            shift, idxs = regex.match(row).groups()
                            shift, idxs = float(shift), [int(idx) for idx in idxs.split(b',')]
                            for atom_idx in idxs:
                                res[mol_idx, atom_idx] = shift
                    mol_idx += 1
        res = pd.Panel(res)
        return res

    def transform(self, inp):
        return super(ChemAxonNMRPredictor, self).transform(inp).T