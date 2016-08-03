#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors

A module concerned with calculating molecular descriptors.
"""

from .physicochemical import PhysicochemicalFeaturizer
from .fingerprints import (AtomPairFeaturizer,
                           MorganFeaturizer, MACCSFeaturizer,
                           TopologicalTorsionFeaturizer, RDKFeaturizer,
                           ErGFeaturizer, ConnectivityInvariantsFeaturizer,
                           FeatureInvariantsFeaturizer)
from .chemaxon import (ChemAxonAtomFeaturizer, ChemAxonFeaturizer, ChemAxonNMRPredictor)
from .atom import (AtomFeaturizer, GraphDistanceTransformer, SpacialDistanceTransformer)

DEFAULTS = {
    'atom': AtomFeaturizer,
    'graph_distance': GraphDistanceTransformer,
    'spacial_distance': SpacialDistanceTransformer,
    'morgan': MorganFeaturizer,
    'atom_pair': AtomPairFeaturizer,
    'topological_torsion': TopologicalTorsionFeaturizer,
    'rdk': RDKFeaturizer,
    'erg': ErGFeaturizer,
    'conn_inv': ConnectivityInvariantsFeaturizer,
    'feat_inv': FeatureInvariantsFeaturizer,
    'physicochemical': PhysicochemicalFeaturizer,
    'chemaxon': ChemAxonFeaturizer,
    'chemaxon_atom': ChemAxonAtomFeaturizer
}

def get(name):
    """ Retrieve a descriptor calculator by name."""

    if isinstance(name, str):
        return DEFAULTS[name]()
    else:
        return name