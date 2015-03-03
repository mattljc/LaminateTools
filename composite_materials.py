class CompositeMaterial():
	# Superclass for all composite material types.
	# The intent is to create a generic material superclass for class testing, along with fundamental functions.

	def __str__(self):
		# Simple string returm
		return self.Name

	def __repr__(self):
		# repr form is the class name with a print of the properties dictionary
		output = self.__class__.__name__ + str(self.__dict__)
		return output

	# add overloads for equals
	# add xml tag output function

class ContinuumMaterial(CompositeMaterial):
	# Defines a 'thick' composite material, with out-of-plane properties. This class is the archetype for all other populated material classes.

	def __init__(self, name=None, E11_in=0, E22_in=0, E33_in=0, Nu12_in=0, Nu13_in=0, Nu23_in=0, G12_in=0, G13_in=0, G23_in=0, ArealDensity_in=0, CPT_in=0):
		# Constructors defaults all numeric properties to zero and strings to null.
		self.Name = name
		self.E11 = E11_in
		self.E22 = E22_in
		self.E33 = E33_in
		self.Nu12 = Nu12_in
		self.Nu13 = Nu13_in
		self.Nu23 = Nu23_in
		self.G12 = G12_in
		self.G13 = G13_in
		self.G23 = G23_in
		self.ArealDensity = ArealDensity_in
		self.CPT = CPT_in

class PlateMaterial(CompositeMaterial):
	# Defines a plate material for use in classic lamination theory. See wiki for where this type is appropriate.

	def __init__(self, name=None, E11_in=0, E22_in=0, Nu12_in=0, G12_in=0, ArealDensity_in=0, CPT_in=0):
		# Constructors defaults all numeric properties to zero and strings to null.
		self.Name = name
		self.E11 = E11_in
		self.E22 = E22_in
		self.Nu12 = Nu12_in
		self.G12 = G12_in
		self.ArealDensity = ArealDensity_in
		self.CPT = CPT_in
