#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.utils.suppress

Class for suppressing C extensions output.
"""

import os

##
# Adapted from:
# http://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
##

class Suppressor(object):
    """ A context manager for doing a "deep suppression" of stdout and stderr.

    It will suppress all print, even if the print originates in a compiled C/Fortran sub-function.

    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    """

    # always have these on the class, so we can have multiple Suppressors with
    # out running out of file descriptors
    null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]

    def __init__(self):
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(Suppressor.null_fds[0], 1)
        os.dup2(Suppressor.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        os.close(self.save_fds[0])
        os.close(self.save_fds[1])

