#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.io.objects

Implementing object serialization and deserialization.
"""

import importlib
import json
import yaml


def read_config(conf):

    """ Deserialize an object from a config dict.

    Args:
        conf (dict):
            The config dict to deseriailize.

    Returns:
        object

    Note:
        `config` is different from `params`, in that it specifies the class.
        The `params` dict is a subdict in `config`.
    """

    (name, params), = conf.items()
    mod, kls = name.rsplit('.', 1)
    m = importlib.import_module(mod)
    return getattr(m, kls).from_params(params)


def write_config(obj):

    """ Serialize an object to a config dict. """

    return obj.to_dict()


def read_json(conf):

    """ Deserialize an object from a json file, filename or str.

    Args:
        json (str or filelike):
            The json file to deserialize.

    Returns:
        object
    """

    if isinstance(conf, str):
        try:
            with open(conf, 'r') as f:
                conf = json.load(f)
        except (FileNotFoundError, OSError):
            conf = json.loads(conf)
    else:
        conf = json.load(conf)
    return read_config(conf)


def write_json(obj, target=None):

    """ Serialize a scikit-chem object as json. """

    return obj.to_json(target)


def read_yaml(conf):

    """ Deserialize an object from a yaml file, filename or str.

    Args:
        yaml (str or filelike):
            The yaml file to deserialize.

    Returns:
        object
    """
    if isinstance(conf, str):
        try:
            with open(conf, 'r') as f:
                conf = yaml.load(f)
        except (FileNotFoundError, OSError):
            conf = yaml.load(conf)

    else:
        conf = yaml.load(conf)
    return read_config(conf)


def write_yaml(obj, target=None):

    """ Serialize a scikit-chem object to yaml. """

    return obj.to_yaml(target)