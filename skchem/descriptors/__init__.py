#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors

A module concerned with calculating molecular descriptors.
"""

from .physicochemical import PhysicochemicalFingerprinter
from .fingerprints import (Fingerprinter, AtomPairFingerprinter,
                           MorganFingerprinter, MACCSKeysFingerprinter,
                           TopologicalTorsionFingerprinter, RDKFingerprinter,
                           ErGFingerprinter, ConnectivityInvariantsFingerprinter,
                           FeatureInvariantsFingerprinter)
from .atom import (AtomFeatureCalculator, GraphDistanceCalculator)

DEFAULTS = {
    'morgan': MorganFingerprinter,
    'atom_pair': AtomPairFingerprinter,
    'topological_torsion': TopologicalTorsionFingerprinter,
    'rdk': RDKFingerprinter,
    'erg': ErGFingerprinter,
    'conn_inv': ConnectivityInvariantsFingerprinter,
    'feat_inv': FeatureInvariantsFingerprinter,
    'physicochemical': PhysicochemicalFingerprinter
}
def get(name):
    """ Retrieve a descriptor calculator by name."""

    if isinstance(name, str):
        return DEFAULTS[name]()
    else:
        return name