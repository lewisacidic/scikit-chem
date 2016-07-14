#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
# skchem.utils.io

IO helper functions for skchem.
"""

import subprocess

def line_count(filename):

    """ Adapted from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python"""
    f = open(filename, 'r')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines

def sdf_count(filename):

    """ Efficiently count molecules in an sdf file.

    Specifically, the function counts the number of times '$$$$\n' occurs in the file.

    This function is quicker than calling subprocesses for sdf files of ~10 000.

    Args:
        filename (str): The filename of the sdf file.

    Returns:
        int: the number of molecules in the file.
    """

    with open(filename, 'rb') as f:
        return sum(1 for l in f if l == b'$$$$\n')