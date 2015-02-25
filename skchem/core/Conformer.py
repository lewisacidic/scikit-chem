import rdkit.Chem
from skchem.core import Point3D

class Conformer(rdkit.Chem.rdchem.Conformer):

    #should use a view, list will do for now
    @property
    def atom_positions(self):
        return [Point3D._from_super(self.GetAtomPosition(i)) for i in range(self.GetNumAtoms())]
    @atom_positions.setter
    def atom_positions(self, value):
        raise NotImplementedError
    
    @property
    def is_3d(self):
        return self.is3D()
    @is_3d.setter
    def is_3d(self, value):
        self.set3D(value)

    def __repr__(self):
        return '<{klass} id="{id}" at {address}>'.format(klass=self.__class__.__name__,
            id=self.GetId(), address=hex(id(self)))

