import rdkit.Geometry.rdGeometry
from skchem.core import ChemicalObject as _ChemicalObject

class Point3D(rdkit.Geometry.rdGeometry.Point3D, _ChemicalObject):

    def to_dict(self, two_d=True):
        if two_d:
            return {"x": round(self.x), "y": round(self.y)}
        else:
            return {"x": round(self.x), "y": round(self.y), "z": round(self.z)}

    def __repr__(self):
        return '<{klass} coords="({x:.2f}, {y:.2f}, {z:.2f})" at {address}>'.format(
            klass=self.__class__.__name__, 
            x=self.x, 
            y=self.y, 
            z=self.z, 
            address=hex(id(self)))

    def __str__(self):
        return '({x:.2f}, {y:.2f}, {z:.2f})'.format(x=self.x, y=self.y, z=self.z)
