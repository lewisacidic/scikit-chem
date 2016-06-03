#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.vis.mol

Module for drawing molecules.
"""

from matplotlib import pyplot as plt
from rdkit.Chem.Draw import DrawingOptions, MolToImage

def draw(mol, quality=1, ax=None):

    """Draw a molecule on a matplotlib axis.

    Args:
        mol (skchem.Mol):
            The molecule to be drawn.
        quality (int):
            The level of quality.  Higher quality takes more time, but will be
            higher quality (so long as matplotlib's savefig.dpi is high enough).

    Returns:
        plt.AxesImage:
            A matplotlib AxesImage object with the molecule drawn.
    """

    if not ax:
        ax = plt.gca()
    ax.grid('off')
    ax.axis('off')

    opts = DrawingOptions()
    opts.dotsPerAngstrom *= quality
    opts.atomLabelFontSize *= quality
    opts.bondLineWidth *= quality

    size = 300 * quality

    img, canvas, drawer = MolToImage(mol, size=(size, size), options=opts, returnCanvas=True)
    canvas.flush()
    return ax.imshow(img, extent=(0, 1, 0, 1))


    return ax
