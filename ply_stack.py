from composite_materials import *

class Laminate(object):
	#
	def __init__(self, plyBook=None, n_count=1, symmetry=False):
		# Test type of plybook, test type of contents
		assert isinstance(plyBook, list)
		for thing in plyBook:
			assert isinstance(thing, Ply)

		# Store symmetry, ncount and apply to ply stack
		self.Symmetry = symmetry
		self.nCount = n_count
		plyBook = plyBook*self.nCount
		if self.Symmetry:
			plyBook += plyBook[::-1] #Fancy slice adding the reverse back to itself
		self.PlyStack = plyBook

	def __str__(self):
		plyNumber = 1
		output = str()
		for ply in self.PlyStack:
			output = output+('Ply #'+str(plyNumber)+' > '+str(ply)+'\n')
			plyNumber += 1
		return output

	def __repr__(self):
		output = '<plybook> \n'
		for ply in self.PlyStack:
			output += ('   '+ply.__repr__()+'\n')
		output += '</plybook>'


class Ply(object):
	#
	def __init__(self, matl=None, orient=0, thk=None):
		# Test material is appropriate type
		assert isinstance(matl, CompositeMaterial)

		self.Material = matl

		# Test ply thickness, use CPT from material if blank
		if thk==None:
			self.Thickness = self.Material.CPT
		else:
			self.Thickness = thk

		self.Orientation = orient

	def __str__(self):
		output = (str(self.Material)+' orient='+str(self.Orientation)+' thickness='+str(self.Thickness))
		return output

	def __repr__(self):
		return str('<ply material=\"'+str(self.Material)+'\" orientation=\"'+str(self.Orientation)+'\" thickness=\"'+str(self.Thickness)+'\" />')
