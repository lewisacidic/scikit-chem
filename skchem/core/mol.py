#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.mol

Defining molecules in scikit-chem.
"""

import rdkit.Chem
import rdkit.Chem.inchi
from rdkit.Chem import AddHs, RemoveHs
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt

import json

from .atom import AtomView
from .bond import BondView
from .conformer import Conformer
from .base import ChemicalObject, PropertyView
from ..utils import Suppressor


class Mol(rdkit.Chem.rdchem.Mol, ChemicalObject):

    """Class representing a Molecule in scikit-chem.

    Mol objects inherit directly from rdkit Mol objects.  Therefore, they
    contain atom and bond information, and may also include properties and
    atom bookmarks.

    Example:
        Constructors are implemented as class methods with the `from_` prefix.

        >>> import skchem
        >>> m = skchem.Mol.from_smiles('CC(=O)Cl'); m # doctest: +ELLIPSIS
        <Mol name="None" formula="C2H3ClO" at ...>

        This is an rdkit Mol:

        >>> from rdkit.Chem import Mol as RDKMol
        >>> isinstance(m, RDKMol)
        True

        A name can be given at initialization:
        >>> m = skchem.Mol.from_smiles('CC(=O)Cl', name='acetyl chloride'); m # doctest: +ELLIPSIS
        <Mol name="acetyl chloride" formula="C2H3ClO" at ...>

        >>> m.name
        'acetyl chloride'

        Serializers are implemented as instance methods with the `to_` prefix.

        >>> m.to_smiles()
        'CC(=O)Cl'

        >>> m.to_inchi()
        'InChI=1S/C2H3ClO/c1-2(3)4/h1H3'

        >>> m.to_inchi_key()
        'WETWJCDKMRHUPV-UHFFFAOYSA-N'

        RDKit properties are accessible through the `props` property:

        >>> m.SetProp('example_key', 'example_value') # set prop with rdkit directly
        >>> m.props['example_key']
        'example_value'

        >>> m.SetIntProp('float_key', 42) # set int prop with rdkit directly
        >>> m.props['float_key']
        42

        They can be set too:

        >>> m.props['example_set'] = 'set_value'
        >>> m.GetProp('example_set') # getting with rdkit directly
        'set_value'

        We can export the properties into a dict or a pandas series:

        >>> m.props.to_series()
        example_key    example_value
        example_set        set_value
        float_key                 42
        dtype: object

        Atoms and bonds are provided in views:

        >>> m.atoms # doctest: +ELLIPSIS
        <AtomView values="['C', 'C', 'O', 'Cl']" at ...>

        >>> m.bonds # doctest: +ELLIPSIS
        <BondView values="['C-C', 'C=O', 'C-Cl']" at ...>

        These are iterable:
        >>> [a.element for a in m.atoms]
        ['C', 'C', 'O', 'Cl']

        The view provides shorthands for some attributes to get these as pandas objects:

        >>> m.atoms.element
        atom_idx
        0     C
        1     C
        2     O
        3    Cl
        dtype: object

        Atom and bond props can also be set:

        >>> m.atoms[0].props['atom_key'] = 'atom_value'
        >>> m.atoms[0].props['atom_key']
        'atom_value'

        The properties for atoms on the whole molecule can be accessed like so:

        >>> m.atoms.props # doctest: +ELLIPSIS
        <MolPropertyView values="{'atom_key': ['atom_value', None, None, None]}" at ...>

        The properties can be exported as a pandas dataframe
        >>> m.atoms.props.to_frame()
                    atom_key
        atom_idx
        0         atom_value
        1               None
        2               None
        3               None
    """

    def __init__(self, *args, **kwargs):

        """
        The default constructor.

        Note:
            This will be rarely used, as it can only create an empty molecule.

        Args:
            *args: Arguments to be passed to the rdkit Mol constructor.
            **kwargs: Arguments to be passed to the rdkit Mol constructor.
        """
        super(Mol, self).__init__(*args, **kwargs)
        self.__two_d = None #set in constructor

    @property
    def name(self):

        """ str: The name of the molecule.

        Raises:
            KeyError"""

        try:
            return self.GetProp('_Name')
        except KeyError:
            return None

    @name.setter
    def name(self, value):

        if value is None:
            self.ClearProp('_Name')
        else:
            self.SetProp('_Name', value)

    @property
    def atoms(self):

        """ List[skchem.Atom]: An iterable over the atoms of the molecule. """

        if not hasattr(self, '_atoms'):
            self._atoms = AtomView(self)
        return self._atoms

    @property
    def bonds(self):

        """ List[skchem.Bond]: An iterable over the bonds of the molecule. """

        if not hasattr(self, '_bonds'):
            self._bonds = BondView(self)
        return self._bonds

    @property
    def mass(self):

        """ float: the mass of the molecule. """

        return CalcExactMolWt(self)

    @property
    def props(self):

        """ PropertyView: A dictionary of the properties of the molecule. """

        if not hasattr(self, '_props'):
            self._props = PropertyView(self)
        return self._props

    @property
    def conformers(self):

        """ List[Conformer]: conformers of the molecule. """

        return [Conformer.from_super(self.GetConformer(i)) \
            for i in range(len(self.GetConformers()))]

    def to_formula(self):

        """ str: the chemical formula of the molecule.

        Raises:
            RuntimeError"""

        # formula may be undefined if atoms are uncertainly typed
        # e.g. if the molecule was initialize through SMARTS
        try:
            return CalcMolFormula(self)
        except RuntimeError:
            raise ValueError('Formula is undefined for {}'.format(self))

    def _two_d(self):

        """ Return a conformer with coordinates in two dimension. """

        if not hasattr(self, '__two_d'):
            self.__two_d = Compute2DCoords(self)
        return self.conformers[self.__two_d]

    def add_hs(self, inplace=False, add_coords=True, explicit_only=False, only_on_atoms=False):
        """

        Args:
            inplace (bool):
                Whether to add Hs to `Mol`, or return a new `Mol`. Default is `False`, return a new `Mol`.
            add_coords (bool):
                Whether to set 3D coordinate for added  Hs. Default is `True`.
            explicit_only (bool):
                Whether to add only explicit Hs, or also implicit ones. Default is `False`.
            only_on_atoms (iterable<bool>):
                An iterable specifying the atoms to add Hs.
        Returns:
            skchem.Mol:
                `Mol` with Hs added.
        """
        if inplace:
            raise NotImplementedError('Inplace addition of Hs is not yet supported.')
        return self.__class__.from_super(AddHs(self, addCoords=add_coords,
                                           onlyOnAtoms=only_on_atoms, explicitOnly=explicit_only))

    def remove_hs(self, inplace=False, sanitize=True, update_explicit=False, implicit_only=False):

        """

        Args:
            inplace (bool):
                Whether to add Hs to `Mol`, or return a new `Mol`. Default is `False`, return a new `Mol`.
            sanitize (bool):
                Whether to sanitize after Hs are removed. Default is `True`.
            update_explicit (bool):
                Whether to update explicit count after the removal. Default is `False`.
            implicit_only (bool):
                Whether to remove explict and implicit Hs, or Hs only. Default is `False`.
        Returns:
            skchem.Mol:
                `Mol` with Hs removed.
        """
        if inplace:
            raise NotImplementedError('Inplace removed of Hs is not yet supported.')
        return self.__class__.from_super(RemoveHs(self, implicitOnly=implicit_only,
                                              updateExplicitCount=update_explicit, sanitize=sanitize))

    def to_dict(self, kind="chemdoodle"):

        """ A dictionary representation of the molecule.

        Args:
            kind (str):
                The type of representation to use.  Only `chemdoodle` is
                currently supported.
                Defaults to 'Chemdoodle'.

        Returns:
            dict:
                dictionary representation of the molecule."""

        if kind == "chemdoodle":
            return self._to_dict_chemdoodle()

        else:
            raise NotImplementedError

    def _to_dict_chemdoodle(self):

        """ Chemdoodle dict representation of the molecule.

        Documentation of the format may be found on the [chemdoodle website](https://web.chemdoodle.com/docs/chemdoodle-json-format/)"""

        atom_positions = [p.to_dict() for p in self._two_d().atom_positions]
        atom_elements = [a.element for a in self.atoms]

        for i, atom_position in enumerate(atom_positions):
            atom_position['l'] = atom_elements[i]

        bonds = [b.to_dict() for b in self.bonds]

        return {"m": [{"a": atom_positions, "b": bonds}]}

    def to_json(self, kind='chemdoodle'):

        """ Serialize a molecule using JSON.

        Args:
            kind (str):
                The type of serialization to use.  Only `chemdoodle` is
                currently supported.

        Returns:
            str: the json string. """

        return json.dumps(self.to_dict(kind=kind))

    def to_inchi_key(self):

        """ The InChI key of the molecule.

        Returns:
            str: the InChI key.

        Raises:
            RuntimeError"""

        if not rdkit.Chem.inchi.INCHI_AVAILABLE:
            raise ImportError("InChI module not available.")

        res = rdkit.Chem.InchiToInchiKey(self.to_inchi())

        if res is None:
            raise RuntimeError("The molecule could not be encoded as InChI key.")

        return res

    def to_binary(self):

        """  Serialize the molecule to binary encoding.

        Args:
            None

        Returns:
            bytes: the molecule in bytes.

        Notes:
            Due to limitations in RDKit, not all data is serialized.  Notably,
            properties are not, so e.g. compound names are not saved."""

        return self.ToBinary()

    @classmethod
    def from_binary(cls, binary):

        """ Decode a molecule from a binary serialization.

        Args:
            binary: The bytes string to decode.

        Returns:
            skchem.Mol: The molecule encoded in the binary."""

        return cls(binary)

    def __repr__(self):
        try:
            formula = self.to_formula()
        except ValueError:
            # if we can't generate the formula, just say it is unknown
            formula = 'unknown'

        return '<{klass} name="{name}" formula="{formula}" at {address}>'.format(
            klass=self.__class__.__name__,
            name=self.name,
            formula=formula,
            address=hex(id(self)))

    def __contains__(self, item):
        if isinstance(item, Mol):
            return self.HasSubstructMatch(item)
        else:
            raise NotImplementedError('No way to check if {} contains {}'.format(self, item))

    def __eq__(self, item):
        if isinstance(item, self.__class__):
            return (self in item) and (item in self)
        else:
            return False

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
        with Suppressor():
            return getattr(rdkit.Chem, 'MolTo' + serializer_name)(self, *args, **kwargs)

    setattr(Mol, 'to_{}'.format(serializer_name).lower() \
        if name_to_bind is None else name_to_bind, serializer)

CONSTRUCTORS = ['Inchi', 'Smiles', 'Mol2Block', 'Mol2File', 'MolBlock', \
                    'MolFile', 'PDBBlock', 'PDBFile', 'Smarts', 'TPLBlock', 'TPLFile']
SERIALIZERS = ['Inchi', 'Smiles', 'MolBlock', 'MolFile', 'PDBBlock', 'Smarts', 'TPLBlock', 'TPLFile']

list(map(bind_constructor, CONSTRUCTORS))
list(map(bind_serializer, SERIALIZERS))
