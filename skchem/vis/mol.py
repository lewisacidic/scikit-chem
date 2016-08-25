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
from matplotlib.widgets import CheckButtons
from mpl_toolkits.mplot3d import Axes3D  # needed to get 3D
import numpy as np

def draw(mol, quality=1, ax=None):

    """Draw a molecule on a matplotlib axis.

    Args:
        mol (skchem.Mol):
            The molecule to be drawn.
        quality (int):
            The level of quality.  Higher quality takes more time, but will
            look better (so long as matplotlib's savefig.dpi is high enough).
        ax (plt.Axes or None):
            An existing axis on which to draw the molecule.

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

    img, canvas, drawer = MolToImage(mol, size=(size, size), options=opts,
                                     returnCanvas=True)
    canvas.flush()
    return ax.imshow(img, extent=(0, 1, 0, 1))


def _plot_sphere(pos=(0, 0, 0), radius=1., ax=None, **kwargs):

    """ Plot a sphere.

    Args:
        pos tuple[float]:
            the centre of the sphere.
        radius (bool):
            The radius of the sphere.
        ax (plt.axes):
            The axes to draw the sphere on (defaults to current axis).

    Note:
        Keyword args are passed to Axes3D.plot_surface.
    """

    if ax is None:
        ax = plt.gca()

    u, v = np.mgrid[0:2 * np.pi:100j], np.mgrid[0:np.pi:100j]

    x = 2 * radius * np.outer(np.cos(u), np.sin(v)) + pos[0]
    y = 2 * radius * np.outer(np.sin(u), np.sin(v)) + pos[1]
    z = 2 * radius * np.outer(np.ones(np.size(u)), np.cos(v)) + pos[2]

    return ax.plot_surface(x, y, z, **kwargs)


_skchem_mpl3d_check = None


def draw_3d(m, conformer_id=-1, label_atoms=None):

    """ Draw a molecule in three dimensions.

    Args:
        conformer_id (int):
            The id of the conformer to draw.
        label_atoms (bool):
            Whether to label the atoms (this can be toggled in interactive
            mode):

    Returns:
        plt.figure

    Note:
        This works great in the notebook with `%matplotlib notebook`.
    """

    # required to stop widgets getting garbage collected
    global _skchem_mpl3d_check

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    pos = m.conformers[conformer_id].positions

    for i, j in m.bonds.atom_idxs:
        ax.plot(*[(pos[i, d], pos[j, d]) for d in range(3)], c='black',
                zorder=1)

    for atom, p in zip(m.atoms, pos):
       _plot_sphere(p, radius=0.25 * atom.covalent_radius, ax=ax,
                    color=atom.hexcode, linewidth=0, zorder=10)
    x, y, z = p + 0.2
    labels = [ax.text(x, y, z, i, visible=False) for i, p in enumerate(pos)]

    # set limits
    lo, hi = min(pos.flatten()) - 1, max(
        pos.flatten()) + 1.5  # 1.5 angstrom outside bording box
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_zlim(lo, hi)
    ax.set_aspect('equal')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    rax = plt.axes([0.15, 0.1, 0.1, 0.15])

    _skchem_mpl3d_check = CheckButtons(rax, ('show labels', 'show axes'),
                                       (False, True))

    def toggle(label):
        if label == 'show labels':
            for l in labels:
                l.set_visible(not l.get_visible())

        if label == 'show axes':
            ax._axis3don = not ax._axis3don  # easiest way

        plt.draw()

    _skchem_mpl3d_check.on_clicked(toggle)

    return ax
