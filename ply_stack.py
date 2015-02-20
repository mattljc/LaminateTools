from composite_materials import *

class Laminate():
	#
	def __init__(self, plyBook=None):
		# Test type of plybook, test type of contents
		if isinstance(plyBook, list):
			for thing in plyBook:
				if not isinstance(thing, Ply):
					raise TypeError('item in plyBook is not type Ply')
		else:
			raise TypeError('plyBook is not type List')

		self.PlyStack = plyBook

	def __str__(self):
		plyNumber = 1
		output = str()
		for ply in self.PlyStack:
			output = output+('Ply #'+str(plyNumber)+' > '+str(ply)+'\n')
			plyNumber += 1
		return output

class Ply():
	#

	def __init__(self, matl=None, orient=0, thk=None):
		# Test material is appropriate type
		if isinstance(matl, RealCompositeMaterial):
			self.Material = matl
		else:
			raise TypeError('material is not a composite-material type')

		# Test ply thickness, use CPT from material if blank
		if thk==None:
			self.Thickness = self.Material.CPT
		else:
			self.Thickness = thk

		self.Orientation = orient

	def __str__(self):
		output = (str(self.Material)+' , '+str(self.Orientation)+' , '+str(self.Thickness))
		return output

# Self test code
if __name__ == '__main__':

	# These statements should all work.
	fakeMaterial = RealCompositeMaterial(name='FakieCF', E11_in=1, E22_in=1, E33_in=1, Nu12_in=1, Nu13_in=1, Nu23_in=1, G12_in=1, G13_in=1, G23_in=1, ArealDensity_in=1, CPT_in=1)
	fakePly = Ply(matl = fakeMaterial, orient=45)
	fakeLaminate = Laminate([fakePly, fakePly, fakePly, fakePly])
	print(fakeLaminate)

	# These statements should not
