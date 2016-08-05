#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.pipeline.pipeline

Module implementing pipelines.
"""

def is_transform_filter(obj):
    """ Whether an object is a TransformFilter (by duck typing). """
    return hasattr(obj, 'transform_filter')

def is_filter(obj):
    """ Whether an object is a Filter (by duck typing). """
    return hasattr(obj, 'filter')

def is_transformer(obj):
    """ Whether an object is a Transformer (by duck typing). """
    return hasattr(obj, 'transform')


class Pipeline(object):

    """ Pipeline object. Applies filters and transformers in sequence. """

    def __init__(self, objects):
        self.objects = objects

    def transform_filter(self, mols, y=None):
        for obj in self.objects:
            if is_transform_filter(obj):
                mols = obj.transform_filter(mols)
            elif is_filter(obj):
                mols = obj.filter(mols)
            elif is_transformer(obj):
                mols = obj.transform(mols)
            else:
                raise NotImplementedError('Cannot apply {}.'.format(obj))
        return mols if y is None else (mols, y.reindex(mols.index))

