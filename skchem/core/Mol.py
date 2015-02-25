import rdkit.Chem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as _molecular_formula
import json
from skchem.core import *
from skchem.core import ChemicalObject as _ChemicalObject

class Mol(rdkit.Chem.rdchem.Mol, _ChemicalObject):

    @property
    def name(self):
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
        return [Atom._from_super(self.GetAtomWithIdx(i)) for i in range(self.GetNumAtoms())]
    @atoms.setter
    def atoms(self, value):
        raise NotImplementedError

    @property
    def bonds(self):
        return [Bond._from_super(self.GetBondWithIdx(i)) for i in xrange(self.GetNumBonds())]
    @bonds.setter
    def bonds(self, value):
        raise NotImplementedError

    @property
    def props(self):
        return {i: self.GetProp(i) for i in self.GetPropNames()}
    @props.setter
    def props(self, value):
        map(self.ClearProp, self.GetPropNames())
        map(lambda v: self.SetProp(v, value[v]), value)

    @property
    def conformers(self):
        return [ Conformer._from_super(self.GetConformer(i)) for i in range(len(self.GetConformers()))]
    @conformers.setter
    def conformers(self, value):
        raise NotImplementedError

    def to_formula(self):
        return _molecular_formula(self)

    def _2D(self):
        if not hasattr(self, '__2D'):
            self.__2D = Compute2DCoords(self)
        return self.conformers[self.__2D]

    def to_dict(self, kind="chemdoodle"):

        if kind == "chemdoodle":
            return self._to_dict_chemdoodle()

        else:
            raise NotImplementedError

    def _to_dict_chemdoodle(self):

        aps = map(lambda p: p.to_dict(), self._2D().atom_positions)
        ats = map(lambda a: a.element, self.atoms)
        for i, ap in enumerate(aps):
            ap['l'] = ats[i]
        bs = map(lambda b: b.to_dict(), self.bonds)
            
        return {"m": [{"a": aps, "b": bs}]}

    def to_json(self):

        return json.dumps(self.to_dict())

    def to_inchi_key(self, *args, **kwargs):

        return rdkit.Chem.InchiToInchiKey(self.to_inchi())

    def __repr__(self):
        return '<{klass} name="{name}" formula="{formula}" at {address}>'.format(
            klass=self.__class__.__name__,
            name=self.name,
            formula=self.to_formula(),
            address=hex(id(self)))

    def _repr_javascript(self):
            return """"""

    def __str__(self):
        return self.to_smiles()

def bind_constructor(constructor_name, to=None):

    @classmethod
    def constructor(self, m, name=None, *args, **kwargs):
        m = Mol(getattr(rdkit.Chem, 'MolFrom' + constructor_name)(m, *args, **kwargs)) 
        m.name = name
        return m     

    setattr(Mol, 'from_{}'.format(constructor_name).lower() if to is None else to, constructor)

def bind_serializer(serializer_name, to=None):

    def serializer(self, *args, **kwargs):
        return getattr(rdkit.Chem, 'MolTo' + serializer_name)(self, *args, **kwargs)

    setattr(Mol, 'to_{}'.format(serializer_name).lower() if to is None else to, serializer)

map(bind_constructor, ['Inchi', 'Smiles', 'Mol2Block', 'Mol2File', 'MolBlock', 'MolFile', 'PDBBlock', 'PDBFile', 'Smarts', 'TPLBlock', 'TPLFile'])
map(bind_serializer, ['Inchi', 'Smiles', 'MolBlock', 'PDBBlock', 'Smarts', 'TPLBlock', 'TPLFile'])

