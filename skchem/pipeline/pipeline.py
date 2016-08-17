#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.pipeline.pipeline

Module implementing pipelines.
"""

from ..utils import yaml_dump, json_dump
from ..io import read_config, read_json, read_yaml


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

    def get_params(self):
        return {'objects': [obj.to_dict() for obj in self.objects]}

    def to_dict(self):
        """ Return a dictionary representation of the object."""

        return {'skchem.pipeline.pipeline.Pipeline': self.get_params()}

    @classmethod
    def from_params(cls, params):
        """ Create a instance from a params dictionary. """
        return cls([read_config(conf) for conf in params['objects']])

    def to_json(self, target=None):

        """ Serialize the object as JSON.

        Args:
            target (str or file-like):
                A file or filepath to serialize the object to.  If `None`,
                return the JSON as a string.

            Returns:
                None or str
        """

        return json_dump(self.to_dict(), target)

    def to_yaml(self, target=None):

        """ Serialize the object as YAML.

        Args:
            target (str or file-like):
                A file or filepath to serialize the object to.  If `None`,
                return the YAML as a string.

            Returns:
                None or str
        """

        return yaml_dump(self.to_dict(), target)

    def copy(self):
        """ Return a copy of this object. """
        return self.__class__(self.get_params())

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

