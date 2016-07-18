#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import re

def camel_to_snail(s):
    return re.sub('(?!^)([A-Z]+)', r'_\1', s).lower()

def free_to_snail(s):
    return s.strip().lower().replace(' ', '_')