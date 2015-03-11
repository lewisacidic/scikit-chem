class ChemicalObject(object):
    
    @classmethod
    def from_super(cls, obj):
        obj.__class__ = cls
        return obj

    pass