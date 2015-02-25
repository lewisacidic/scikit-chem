import rdkit.Chem
from skchem.core import Atom
from skchem.core import ChemicalObject as _ChemicalObject

class Bond(rdkit.Chem.rdchem.Bond, _ChemicalObject):

    @property
    def order(self):
        return self.GetBondTypeAsDouble()
    @order.setter
    def order(self, value):
        raise NotImplementedError

    @property
    def atoms(self):
        return [Atom._from_super(self.GetBeginAtom()), Atom._from_super(self.GetEndAtom())]
    @atoms.setter
    def atoms(self, value):
        raise NotImplementedError
    
    def draw(self):
        return '{}{}{}'.format(self.atoms[0].element, '-' if self.order == 1 else self.GetSmarts(), self.atoms[0].element)
    
    def to_dict(self):
        return {"b": self.GetBeginAtomIdx(), "e": self.GetEndAtomIdx(), "o": self.order}
        
    def __repr__(self):
        return '<{klass} type="{bond}" at {address}>'.format(klass=self.__class__.__name__, 
            bond=self.draw(),
            address=hex(id(self)))

    def __str__(self):
        return self.draw()