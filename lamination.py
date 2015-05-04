import materials
from copy import deepcopy

class Laminate(object):
	# A laminate is an ordered list of plies. The order is assumed to
	# start from the tool side as is conventional, although this has
	# little effect on the operation of the code.

	def __init__(self, plyBook=None, n_count=1, symmetry=False):
		# Test type of plybook, test type of contents
		assert isinstance(plyBook, list)
		for thing in plyBook:
			assert isinstance(thing, Ply)

		# Store symmetry, ncount and apply to ply stack.
		# Do not hard code symmetry unless absolutely necessary. Other
		# functions will check this boolean to see if symmetry-only
		# tricks can be used.
		self.Symmetry = symmetry
		self.nCount = n_count

		for ct in range(self.nCount-1):
			temp_list = list()
			for ply in plyBook:
				temp_ply = deepcopy(ply)
				temp_list += [temp_ply]
			plyBook = plyBook + temp_list

		#Fancy slice adding the reverse back to itself
		if self.Symmetry:
			temp_list = list()
			for ply in plyBook[::-1]:
				temp_ply = deepcopy(ply)
				temp_list += [temp_ply]
			plyBook += temp_list

		self.PlyStack = plyBook

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
		pass


class Ply(object):
	# A ply is a single layer of a laminated composite material.
	# Orientation is measured between the logical ply 1-direction and
	# the global x-direction.
	def __init__(self, definitionDict):
		# Test material is appropriate type
		assert isinstance(definitionDict['matl'], materials.Materials)
		try:
			self.Material = definitionDict['matl']
			self.Thickness = definitionDict['thk']
			self.Orientation = definitionDict['orient']
		except KeyError:
			raise KeyError('Check properties input, minimum information not provided')

	def __str__(self):
		# The string output is functionally the same as the representation
		return self.__repr__()

	def __repr__(self):
		output  = '--Ply--\n'
		output += 'material = '+self.Material.__repr__()+'\n'
		output += 'orientation = '+str(self.Orientation)+'\n'
		output += 'thickness = '+str(self.Thickness)
		return output

	def toXML(self):
		pass
