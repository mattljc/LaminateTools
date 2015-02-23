class RealCompositeMaterial():
	# Defines a real composite material. Defaults all properties to zero

	def __init__(self, name=None, E11_in=0, E22_in=0, E33_in=0, Nu12_in=0, Nu13_in=0, Nu23_in=0, G12_in=0, G13_in=0, G23_in=0, ArealDensity_in=0, CPT_in=0):
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

	def __str__(self):
		return self.Name

	def __repr__(self):
		return str('<materialProps name=\"',self.Name,'\" E11=\"',str(self.E11),'\" E22=\"',str(self.E22),'\" E33=\"',str(self.E33),
		'\" nu12=\"',str(self.Nu12),'\" nu13=\"',str(self.Nu13),'\" nu23=\"',str(self.Nu23),
		'\" G12=\"',str(self.G12),'\" G13=\"',str(self.G13),'\" G23=\"',str(self.G23),
		'\" CPT=\"',str(self.CPT),'\" arealDensity=\"',str(self.ArealDensity),'\" />')

# Save these two for a later day

class UnidirectionalFabric(RealCompositeMaterial):

	def __init__(self):
		self.Name = 'name'

class PlainWeaveFabric(RealCompositeMaterial):

	def __init__(self):
		self.Name = 'name'