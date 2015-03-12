import rdkit.Chem as _Chem
from skchem.core import ChemicalObject as _ChemicalObject

class Atom(_Chem.rdchem.Atom, _ChemicalObject):

    @property
    def element(self):
        return self.GetSymbol()
    @element.setter
    def element(self, value):
        raise NotImplementedError

    @property
    def props(self):
        return {i: self.GetProp() for i in self.GetProps()}
    @props.setter
    def props(self, value):
        map(self.ClearProp, self.GetPropNames())
        map(lambda k: self.SetProp(k, value[k]), value)

    def __repr__(self):
        return '<{klass} element="{element}" at {address}>'.format(
            klass=self.__class__.__name__,
            element=self.element,
            address=hex(id(self))
            )

    def __str__(self):
        return self.element
