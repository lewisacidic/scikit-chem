#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
# skchem.utils.progress

Module implementing progress bars.
"""

import progressbar


class NamedProgressBar(progressbar.ProgressBar):

    def __init__(self, name=None, **kwargs):
        kwargs.setdefault('redirect_stderr', True)
        super(NamedProgressBar, self).__init__(**kwargs)
        self.name = name

    def default_widgets(self):
        return [self.name, ': '] + \
               super(NamedProgressBar, self).default_widgets()


class DummyProgressBar(object):

    def __init__(self, *args, **kwargs):
        pass

    def update(self, val):
        pass

    def finish(self):
        pass

    def __call__(self, obj):
        return obj

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass
