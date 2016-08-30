#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.features.descriptors.decorators

Decorators for descriptors in scikit-chem.
"""
import inspect
from functools import wraps
from collections import OrderedDict, defaultdict

import pandas as pd


def requires_h_depleted(func):

    """ Decorate a function that requires an h-depleted graph.

    This will check if the molecule argument is h-depleted, and will
    memoize a depleted version if it is not, and pass it to the func.


    """

    @wraps(func)
    def inner(mol, *args, **kwargs):

        # if is already h depleted
        if (mol.atoms.atomic_number == 1).sum() == 0:
            return func(mol, *args, **kwargs)

        if not hasattr(mol, '_h_depleted'):
            mol._h_depleted = mol.remove_hs()

        return func(mol._h_depleted, *args, **kwargs)

    return inner


def requires_h_filled(func):

    """ Decorate a function that requires a h-filled graph.

    This will check if the molecule argument is h-filled, and will
    memoize a filled version if it is not, and pass it to the func.

     Note:
        This decorator should be used first if in combination with dMat etc.

     """

    @wraps(func)
    def inner(mol, *args, **kwargs):
        # if is already h enriched
        if mol.atoms.n_total_hs.sum() == 0:
            return func(mol, *args, **kwargs)

        # if not, memoize the enriched one and pass it
        if not hasattr(mol, '_h_enriched'):
            mol._h_enriched = mol.add_hs()

        return func(mol._h_enriched, *args, **kwargs)

    return inner


class Cache(object):

    """ Function cache."""

    def __init__(self):
        self.cached = {}

    @staticmethod
    def extract_kwargs(func):

        # py2 compat
        if pd.compat.PY2:
            spec = inspect.getargspec(func)
            if spec.defaults:
                kwds = spec.args[-len(spec.defaults):]
                kwds = OrderedDict(zip(kwds, spec.defaults))
            else:
                kwds = OrderedDict()

        else:

            kwds = OrderedDict((k, v.default) for k, v in
                               inspect.signature(func).parameters.items()
                               if v.default != inspect._empty)

        return kwds

    def __call__(self, func):

        """ Create a decorator to identify a function as returning a cached
        value.

        This can be used for objects  that are nontrivial to generate, but
        are used by many functions.
        """

        name = func.__name__

        # get the key word arguments and the default values of the function

        kwds = self.extract_kwargs(func)

        @wraps(func)
        def inner(mol, *args, **kwargs):

            # py2 compat
            force = kwargs.pop('force', False)

            self.setup_cache(mol)

            # get the full set of keywords to use, including defaults
            kwds.update(kwargs)
            kw_to_save = tuple(sorted(kwds.items()))

            # call function if it hasn't already been called
            # with required arguments, or if told to.
            if force or name not in self.cached.keys() or \
                    kw_to_save not in mol.cache.get(name, {}).keys():

                res = func(mol, *args, **kwargs)

                # cache the value with the args used.
                mol.cache[name].update({kw_to_save: res})

            # return the cached value
            return mol.cache[name][kw_to_save]

        self.cached[name] = inner, tuple(kwds.keys())

        return inner

    def inject(self, *args_to_inject):

        """ Create a decorator that will inject cached values as arguments.

        Args:
            args (list<str>):
                A list of cachable requirements for this function.
        """

        def outer(func):

            # extract the defaults for the func
            kwds = self.extract_kwargs(func)

            @wraps(func)
            def inner(mol, *args, **kwargs):

                # augment with the keywords from the function
                kwds.update(kwargs)
                self.setup_cache(mol)

                # look up cached values, or produce them if not.
                # inject the cached values

                args_supp = ()

                for arg in args_to_inject:
                    # lookup function to inject
                    inj_func, params = self.cached[arg.__name__]

                    # get the kwargs required
                    inj_kwargs = {param: kwds[param] for param in params
                                  if param in kwds.keys()}

                    # get a hashable representation of the kwargs
                    immut = tuple(sorted(inj_kwargs.items()))

                    # retrieve the cached result (or None if not yet cached)
                    res = mol.cache.get(arg.__name__, {}).get(immut, None)

                    # calculate and cache result
                    if res is None:
                        res = inj_func(mol, **inj_kwargs)

                    # add to injected args
                    args_supp += (res,)

                # put injected args at start of arg list
                args = args_supp + args

                return func(mol, *args, **kwargs)

            return inner

        return outer

    @staticmethod
    def setup_cache(mol):

        """ Set up a cache on e.g. a `Mol`. """

        if not hasattr(mol, 'cache'):
            mol.cache = defaultdict(dict)

    @staticmethod
    def teardown_cache(mol):

        """ Tear down a cache on e.g. a `Mol`. """

        if hasattr(mol, 'cache'):
            del mol.cache

cache = Cache()
