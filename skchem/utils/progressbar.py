#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
# skchem.utils.progressbar

Module implementing progressbars.
"""

import progressbar

class NamedProgressBar(progressbar.ProgressBar):

    def __init__(self, name=None, **kwargs):
        super(NamedProgressBar, self).__init__(**kwargs)
        self.name = name

    def default_widgets(self):
        return [self.name, ': '] + super(NamedProgressBar, self).default_widgets()
