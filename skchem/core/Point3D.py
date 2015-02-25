import rdkit.Geometry.rdGeometry

class Point3D(rdkit.Geometry.rdGeometry.Point3D):

    @classmethod
    def _from_super(self, sup):
        sup.__class__ = Point3D
        return sup

    def to_dict(self, two_d=True):
        if two_d:
            return {"x": round(self.x), "y": round(self.y)}
        else:
            return {"x": round(self.x), "y": round(self.y), "z": round(self.z)}

    def __repr__(self):
        return '<{klass} coords="({x:.2f}, {y:.2f}, {z:.2f})">'.format(klass=self.__class__.__name__, x=self.x, y=self.y, z=self.z)
