#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.vis.atom

Module for atom contribution visualization.
"""

from rdkit.Chem.Draw import MolToImage, DrawingOptions

import numpy as np
from matplotlib import pyplot as plt


def plot_weights(mol, weights, quality=1, l=0.4, step=50, levels=20, contour_opacity=0.5, cmap='RdBu', ax=None, **kwargs):
    """ Plot weights as a sum of gaussians across a structure image.

    Args:
        mol (skchem.Mol):
            Molecule to visualize weights for.
        weights (iterable<float>):
            Array of weights in atom index order.
        l (float):
            Lengthscale of gaussians to visualize as a multiple of bond length.
        steps (int):
            Size of grid edge to calculate the gaussians.
        levels (int):
            Number of contours to plot.
        contour_opacity (float):
            Alpha applied to the contour layer.
        ax (plt.axis):
            Axis to apply the plot to. Defaults to current axis.
        cmap (plt.cm):
            Colormap to use for the contour.
        **kwargs:
            Passed to contourf function.

    Returns:
        matplotlib.AxesSubplot: The plot.
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
    coords = np.array([[i / size, 1 - j / size] for k, (i, j) in list(drawer.atomPs.values())[0].items()])

    b = mol.bonds[0]
    begin, end = b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx()
    length = np.linalg.norm(coords[end] - coords[begin])

    x = np.linspace(0, 1, 500)
    y = np.linspace(0, 1, 500)
    x, y = np.meshgrid(x, y)

    def gaussian(x, y, mu=np.zeros(2), sigma=np.identity(2), size=50):
        return (1 / (2 * np.pi * sigma[0, 0] * sigma[1, 1]) * np.exp(-((x - mu[0]) ** 2 / (2 * sigma[0, 0] ** 2)
         + (y - mu[1]) ** 2 / (2 * sigma[1, 1] ** 2))))

    if not np.max(weights) == np.min(weights) == 0:
        z = sum([w * gaussian(x, y, mu, sigma=l * length * np.identity(2)) for mu, w in zip(coords, weights)])
        v = np.max((np.abs(z.min()), np.abs(z.max())))
    else:
        z = np.zeros(x.shape)
        v = 1

    if z.min() >= 0:
        levels = int(levels/2)
    cf = ax.contourf(x, y, z, levels, alpha=contour_opacity, extent=(0, 1, 0, 1), vmin=-v, vmax=v, cmap=cmap, **kwargs)

    ax.imshow(img, extent=(0, 1, 0, 1))
    return ax
