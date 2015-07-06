from copy import deepcopy
import property_interface

"""Defines the two fundamental types required for all laminate analysis.
"""

class Ply(object):
    """A ply is a single layer of material.

    A ply has a material (must be a material object), thickness and
    orientation. The orientation is measured relative to the logical
    x-direction (convention holds this as the longest dimension) of the part.
    """

    def __init__(self, input_dict):
        assert isinstance(input_dict['matl'], property_interface.Material)
        try:
            self.Material = input_dict['matl']
            self.Thickness = float(input_dict['thk'])
            self.Orientation = float(input_dict['orient'])
        except KeyError:
            raise KeyError('Check ply input, minimum information not provided')

    def __str__(self):
        # The string output is pretty much the same as the representation.
        return self.__repr__()

    def __repr__(self):
        output  = '--Ply--\n'
        output += 'material = '+self.Material.__repr__()+'\n'
        output += 'orientation = '+str(self.Orientation)+'\n'
        output += 'thickness = '+str(self.Thickness)
        return output

    def toXML(self):
        raise NotImplementedError


class Laminate(object):
    """A laminate is an ordered sequence of plies.

    A laminate has a list of ply objects and stores its symmetry state. This
    class also provides functions to build laminates that are symetric and/or
    with repeating units. It is assumed that the list starts with the
    tool side.
    """

    def __init__(self, ply_book=None, n_count=1, symmetry=False):
        assert isinstance(ply_book, list)
        for thing in ply_book:
            assert isinstance(thing, Ply)

        # Do not hard code symmetry unless absolutely necessary. Other
        # functions will check this boolean to see if symmetric-only
        # tricks can be used.
        self.Symmetry = bool(symmetry)
        self.nCount = int(n_count)

        # Can't use list operators because deep copies are needed.
        for ct in range(self.nCount-1):
            temp_list = list()
            for ply in ply_book:
                temp_ply = deepcopy(ply)
                temp_list += [temp_ply]
                ply_book = ply_book + temp_list

        # Fancy slice adding the reverse back to itself
        if self.Symmetry:
            temp_list = list()
            for ply in ply_book[::-1]:
                temp_ply = deepcopy(ply)
                temp_list += [temp_ply]
            ply_book += temp_list

        self.PlyStack = ply_book

    def __str__(self):
        # The string output is functionally the same as the representation
        return self.__repr__()

    def __repr__(self):
        output  = '--Laminate--\n'
        output += 'symmetry? '+str(self.Symmetry)+'\n'
        for ply in self.PlyStack:
            output += ply.__repr__()+'\n'
        return output

    def toXML(self):
        raise NotImplementedError
