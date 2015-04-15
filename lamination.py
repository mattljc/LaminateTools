import materials

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
		plyBook = plyBook*self.nCount

		#Fancy slice adding the reverse back to itself
		if self.Symmetry:
			plyBook += plyBook[::-1]

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
	def __init__(self, matl=None, orient=0, thk=None):
		# Test material is appropriate type
		#assert isinstance(matl, materials.Materials)
		self.Material = matl
		# Test ply thickness, use CPT from material if blank
		if thk==None:
			self.Thickness = self.Material.CPT
		else:
			self.Thickness = thk
		self.Orientation = orient

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
