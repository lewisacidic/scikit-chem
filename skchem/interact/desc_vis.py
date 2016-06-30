#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

import matplotlib.pyplot as plt

from .. import descriptors
from .. import core
from .. import vis
from ipywidgets import Dropdown, Text, VBox, HBox, Valid, HTML
from IPython import get_ipython
from IPython.display import clear_output, display


class Visualizer(object):
    def __init__(self, fper='morgan', smiles='c1ccccc1O', dpi=200):

        self.initialize_ipython()

        if isinstance(fper, str):
            self.fper = descriptors.get(fper)
        else:
            self.fper = fper

        self.smiles_input = Text(smiles, description='smiles')
        self.smiles_input.on_submit(self.update_smiles)
        self.smiles_input.observe(self.typing)

        self.valid = Valid(True)

        self.dropdown = Dropdown(options=[], description='bit')
        self.dropdown.observe(self.plot)

        self.dpi_input = Text(str(dpi), description='dpi')
        self.dpi_input.on_submit(self.plot)

        self.ui = VBox([
            HTML('<h2>Visualizer</h2>'),
            HBox([self.smiles_input, self.valid]),
            self.dropdown,
            self.dpi_input])

        self.update_smiles(None)
        self.display()

    def initialize_ipython(self):
        ipython = get_ipython()
        try:
            ipython.magic('matplotlib inline')
        except:
            pass

    def typing(self, _):
        self.valid.visible = False

    @property
    def dpi(self):
        try:
            return int(self.dpi_input.value)
        except:
            return 50

    @dpi.setter
    def dpi(self, value):
        self.dpi_input.value = str(value)

    def display(self):
        display(self.ui)

    def update_smiles(self, _):
        try:
            self._mol = core.Mol.from_smiles(self.smiles_input.value)
            self.valid.value = True
        except ValueError:
            self.valid.value = False
            return
        finally:
            self.valid.visible = True
        return self.calculate()

    def calculate(self):
        fp = self.fper.transform(self.mol)
        self.fp = fp[fp == 1].index
        self.fpg = self.fper.grad(self.mol).ix[self.fp]
        return self.update_dropdown()

    def update_dropdown(self):
        self.dropdown.options.append(self.fp[0])
        self.dropdown.value = self.fp[0]
        self.dropdown.options = self.fp.tolist()
        return self.plot(self.dropdown.value)

    @property
    def mol(self):
        return self._mol

    @mol.setter
    def mol(self, mol):
        self._mol = mol
        self.smiles_input.value = mol.to_smiles()
        self.calculate()

    @property
    def current_smiles(self):
        return self.smiles_input.value

    @property
    def current_bit(self):
        return self.dropdown.value

    def plot(self, _):
        clear_output()
        plt.clf()
        plt.rcParams['savefig.dpi'] = self.dpi
        vis.plot_weights(self.mol, self.fpg.ix[self.current_bit], quality=4, ax=plt.gca())
