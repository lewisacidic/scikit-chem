#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" A cheminformatics library to integrate with the Scientific Python Stack """

import logging

from . import core
from . import data
from . import descriptors
from . import io
from . import vis
from . import target_prediction
from . import cross_validation
from . import standardizers
from . import interact

from .core import Mol
from .descriptors import MorganFingerprinter
from .io import read_sdf, read_smiles

__version__ = '0.0.5'

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())
