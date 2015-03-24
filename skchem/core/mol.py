#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.core.mol

Defining molecules in scikit-chem.
"""

import rdkit.Chem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

import json

from skchem.core import ChemicalObject, Atom, Bond, Conformer

class Mol(rdkit.Chem.rdchem.Mol, ChemicalObject):

    """
    Class representing a Molecule in scikit-chem.
    """

    def __init__(self, *args, **kwargs):
        super(Mol, self).__init__(*args, **kwargs)
        self.__two_d = None #set in constructor

    @property
    def name(self):

        """ Return the name of the molecule. """

        try:
            return self.GetProp('_Name')
        except KeyError:
            return None

    @name.setter
    def name(self, value):

        """ Set the name of the molecule. """

        if value is None:
            self.ClearProp('_Name')
        else:
            self.SetProp('_Name', value)

    @property
    def atoms(self):

        """ return an iterable over the atoms of the molecule. """

        return [Atom.from_super(self.GetAtomWithIdx(i)) for i in range(self.GetNumAtoms())]

    @atoms.setter
    def atoms(self, value):

        """ Set the atoms of a molecule.  Not implemented. """

        raise NotImplementedError

    @property
    def bonds(self):

        """ Return an iterable over the bonds of the molecule. """

        return [Bond.from_super(self.GetBondWithIdx(i)) for i in xrange(self.GetNumBonds())]

    @bonds.setter
    def bonds(self, value):

        """ Set the bonds of the molecule. Not implemented. """

        raise NotImplementedError

    # use a view to easily set properties?
    @property
    def props(self):

        """ Return a dictionary of the properties of the molecule. """

        return {i: self.GetProp(i) for i in self.GetPropNames()}

    @props.setter
    def props(self, value):

        """ Set the properties of the molecule. Not implemented. """

        raise NotImplementedError

    @property
    def conformers(self):

        """ Return a list of conformers of the molecule. """

        return [Conformer.from_super(self.GetConformer(i)) \
            for i in range(len(self.GetConformers()))]

    @conformers.setter
    def conformers(self, value):

        """ Set the conformers of the molecule.  Not implemented. """

        raise NotImplementedError

    def to_formula(self):

        """ Return the chemical formula of the molecule. """

        return CalcMolFormula(self)

    def _two_d(self):

        """ Return a conformer with coordinates in two dimension. """

        if not hasattr(self, '__two_d'):
            self.__two_d = Compute2DCoords(self)
        return self.conformers[self.__two_d]

    def to_dict(self, kind="chemdoodle"):

        """
        Return a dictionary representation of the molecule.

        TODO

        """

        if kind == "chemdoodle":
            return self._to_dict_chemdoodle()

        else:
            raise NotImplementedError

    def _to_dict_chemdoodle(self):

        """ Return a chemdoodle dict representation of the molecule. """

        atom_positions = [p.to_dict() for p in self._2D().atom_positions]
        atom_elements = [a.element for a in self.atoms]

        for i, atom_position in enumerate(atom_positions):
            atom_position['l'] = atom_elements[i]

        bonds = [b.to_dict() for b in self.bonds]

        return {"m": [{"a": atom_positions, "b": bonds}]}

    def to_json(self):

        """ Return a JSON representation of the molecule. """

        return json.dumps(self.to_dict())

    def to_inchi_key(self, *args, **kwargs):

        """ Return the INCHI key of the molecule. """

        return rdkit.Chem.InchiToInchiKey(self.to_inchi(), *args, **kwargs)

    def __repr__(self):
        return '<{klass} name="{name}" formula="{formula}" at {address}>'.format(
            klass=self.__class__.__name__,
            name=self.name,
            formula=self.to_formula(),
            address=hex(id(self)))

    def _repr_javascript(self):

        """ Rich printing in javascript. """

        return self.to_json()

    def __str__(self):
        return '<Mol: {}>'.format(self.to_smiles())

def bind_constructor(constructor_name, name_to_bind=None):

    """ Bind an (rdkit) constructor to the class """

    @classmethod
    def constructor(_, in_arg, name=None, *args, **kwargs):

        """ The constructor to be bound. """

        m = getattr(rdkit.Chem, 'MolFrom' + constructor_name)(in_arg, *args, **kwargs)
        if m is None:
            raise ValueError('Failed to parse molecule, {}'.format(in_arg))
        m = Mol.from_super(m)
        m.name = name
        return m

    setattr(Mol, 'from_{}'.format(constructor_name).lower() \
        if name_to_bind is None else name_to_bind, constructor)

def bind_serializer(serializer_name, name_to_bind=None):

    """ Bind an (rdkit) serializer to the class """

    def serializer(self, *args, **kwargs):

        """ The serializer to be bound. """
        return getattr(rdkit.Chem, 'MolTo' + serializer_name)(self, *args, **kwargs)

    setattr(Mol, 'to_{}'.format(serializer_name).lower() \
        if name_to_bind is None else name_to_bind, serializer)

CONSTRUCTORS = ['Inchi', 'Smiles', 'Mol2Block', 'Mol2File', 'MolBlock', \
                    'MolFile', 'PDBBlock', 'PDBFile', 'Smarts', 'TPLBlock', 'TPLFile']
SERIALIZERS = ['Inchi', 'Smiles', 'MolBlock', 'PDBBlock', 'Smarts', 'TPLBlock', 'TPLFile']

list(map(bind_constructor, CONSTRUCTORS))
list(map(bind_serializer, SERIALIZERS))

