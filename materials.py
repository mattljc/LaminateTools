import numpy as np

class Materials(object):
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

class Continuum(Materials):
	# Defines a 'thick' composite material, with out-of-plane properties. This class is the archetype for all other populated material classes.

	def __init__(self, propsDict = None):
		# Collect everything from a dictionary. This first batch are required
		try:
			self.Name = propsDict['name']
			self.E11 = propsDict['E11']
			self.E22 = propsDict['E22']
			self.E33 = propsDict['E33']
			self.Nu12 = propsDict['Nu12']
			self.Nu13 = propsDict['Nu13']
			self.Nu23 = propsDict['Nu23']
			self.G12 = propsDict['G12']
			self.G13 = propsDict['G13']
			self.G23 = propsDict['G23']
			self.Density = propsDict['Dens']
			self.CPT = propsDict['CPT']
		except KeyError:
			throw KeyError('Check properties input, minimum information not provided')

		# These properties are not required for basic functionality.
		# CTE
		try:
			self.a1 = propsDict['a1']
			self.a2 = propsDict['a2']
			self.a3 = propsDict['a3']
		except KeyError:
			raise UserWarning('No CTE included, setting all CTE to zero')
			self.a1 = 0
			self.a2 = 0
			self.a3 = 0
		# Dynamics to be implemented
		# CME to be implemented

		# Make Compliance
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

class Plate(Materials):
	# Defines a plate material for use in classic lamination theory. See wiki for where this type is appropriate.

	def __init__(self, propsDict = None):
		# Collect everything from a dictionary. This first batch are required
		try:
			self.Name = propsDict['name']
			self.E11 = propsDict['E11']
			self.E22 = propsDict['E22']
			self.Nu12 = propsDict['Nu12']
			self.G12 = propsDict['G12']
			self.Density = propsDict['dens']
			self.CPT = propsDict['CPT']
		except KeyError:
			throw KeyError('Check properties input, minimum information not provided')

		# These properties are not required for basic functionality.
		# CTE
		try:
			self.a1 = propsDict['a1']
			self.a2 = propsDict['a2']
		except KeyError:
			raise UserWarning('No CTE included, setting all CTE to zero')
			self.a1 = 0
			self.a2 = 0

		# Make compliance
		s11 = 1/self.E11
		s12 = -self.Nu12/self.E11
		s22 = 1/self.E22
		s66 = 1/self.G12
		self.Compliance = np.matrix([ \
		[s11, s12, 0  ], \
		[s12, s22, 0  ], \
		[0  , 0  , s66]])

		# Make invariants
		Q = self.Compliance.I
		self.U1 = (Q[0,0] + Q[1,1])*3/8 + Q[0,1]/4 + Q[2,2]/2
		self.U2 = (Q[0,0] - Q[1,1])/2
		self.U3 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 - Q[2,2]/2
		self.U4 = (Q[0,0] + Q[1,1])/8 + Q[0,1]*3/4 - Q[2,2]/2
		self.U5 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 + Q[2,2]/2

	def __repr__(self):
		# Overiding until I can figure out a nice generalized way to do this -mljc
		output = ('E11={e11:.3e} E22={e22:.3e} Nu12={nu12:.3g} G12={g12:.3e} ArealDensity={dens:.3g} CPT={cpt:.3g}').format(\
		e11=self.E11,e22=self.E22,nu12=self.Nu12,g12=self.G12,dens=self.ArealDensity,cpt=self.CPT)
		return output

class Beam(Materials):

	def __init__(self):
		raise NotImplementedError

class Isotropic(Materials):

	def __init__(self):
		raise NotImplementedError
