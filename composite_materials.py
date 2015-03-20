import numpy as np

class CompositeMaterial(object):
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

	def __init__(self, name=None, E11_in=0, E22_in=0, E33_in=0, Nu12_in=0, Nu13_in=0, Nu23_in=0, G12_in=0, G13_in=0, G23_in=0, a1_in=0, a2_in=0, a3_in=0, ArealDensity_in=0, CPT_in=0):
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
		self.a1 = a1_in
		self.a2 = a2_in
		self.a3 = a3_in
		self.ArealDensity = ArealDensity_in
		self.CPT = CPT_in

		s11 = 1/self.E11
		s22 = 1/self.E22
		s33 = 1/self.E33
		s12 = -self.Nu12/self.E11
		s13 = -self.Nu13/self.E11
		s23 = -self.Nu23/self.E22
		s66 = 1/(2*self.G12)
		s55 = 1/(2*self.G13)
		s44 = 1/(2*self.G23)

		self.Compliance = np.matrix([ \
		[s11, s12, s13, 0  , 0  , 0  ], \
		[s12, s22, s23, 0  , 0  , 0  ], \
		[s13, s23, s33, 0  , 0  , 0  ], \
		[0  , 0  , 0  , s44, 0  , 0  ], \
		[0  , 0  , 0  , 0  , s55, 0  ], \
		[0  , 0  , 0  , 0  , 0  , s66]])

	def __repr__(self):
		# Overiding until I can figure out a nice generalized way to do this -mljc
		output = ('E11={e11:.3e} E22={e22:.3e} E33={e33:.3e} Nu12={nu12:.3g} Nu13={nu13:.3g} Nu23={nu23:.3g} G12={g12:.3e} G13={g13:.3e} G23={g23:.3e} ArealDensity={dens:.3g} CPT={cpt:.3g}').format(\
		e11=self.E11,e22=self.E22,e33=self.E33,nu12=self.Nu12,nu13=self.Nu13,nu23=self.Nu23,g12=self.G12,g13=self.G13,g23=self.G23,dens=self.ArealDensity,cpt=self.CPT)
		return output

class PlateMaterial(CompositeMaterial):
	# Defines a plate material for use in classic lamination theory. See wiki for where this type is appropriate.

	def __init__(self, name=None, E11_in=0, E22_in=0, Nu12_in=0, G12_in=0, a1_in=0, a2_in=0, ArealDensity_in=0, CPT_in=0):
		# Constructors defaults all numeric properties to zero and strings to null.
		self.Name = name
		self.E11 = E11_in
		self.E22 = E22_in
		self.Nu12 = Nu12_in
		self.G12 = G12_in
		self.a1 = a1_in
		self.a2 = a2_in
		self.ArealDensity = ArealDensity_in
		self.CPT = CPT_in

		s11 = 1/self.E11
		s12 = -self.Nu12/self.E11
		s22 = 1/self.E22
		s66 = 1/self.G12

		self.Compliance = np.matrix([ \
		[s11, s12, 0  ], \
		[s12, s22, 0  ], \
		[0  , 0  , s66]])

		Q = self.Compliance.I
		self.U1 = 3*(Q[0,0] + Q[1,1])/8 + Q[0,1]/4 + Q[2,2]/2
		self.U2 = (Q[0,0] - Q[1,2])/2
		self.U3 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 - Q[2,2]/2
		self.U4 = (Q[0,0] + Q[1,1])/8 + 3*Q[0,1]/4 - Q[2,2]/2
		self.U5 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 + Q[2,2]/2

	def __repr__(self):
		# Overiding until I can figure out a nice generalized way to do this -mljc
		output = ('E11={e11:.3e} E22={e22:.3e} Nu12={nu12:.3g} G12={g12:.3e} ArealDensity={dens:.3g} CPT={cpt:.3g}').format(\
		e11=self.E11,e22=self.E22,nu12=self.Nu12,g12=self.G12,dens=self.ArealDensity,cpt=self.CPT)
		return output
