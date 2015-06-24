import numpy as np
import warnings

class Materials(object):
	# Superclass for all composite material types.
	# The intent is to create a generic material superclass for class
	# testing, along with fundamental functions.

	def __init__(self, propsDict = None):
		assert isinstance(propsDict, dict)
		self.InputDict = propsDict

	def __str__(self):
		# Simple name string returm
		return self.Name

	def __repr__(self):
		# repr form is the class name with a print out of the properties
		output = self.__class__.__name__ +'  ' + str(self.InputDict)
		return output

	# Implementation of logical operators is meaningless beyond equals
	# and not equals. For safety, the user should only use 'is' operator
	# instead of logicals which may have unintended side effects.
	#
	# At this stage do not attempt to mix material types, even thought
	# this is logical in some analysis types.

class Isotropic(Materials):
	# A classic isotropic material

	def __init__(self, propsDict = None):
		Materials.__init__(propsDict)
		self.Name = self.InputDict['name']
		self.E = self.InputDict['E']
		self.Nu = self.InputDict['Nu']
		self.G = self.E / (2*(1+self.Nu))
		self.Density = self.InputDict['Dens']

class Continuum(Materials):
	# Defines a 'thick' composite material, with out-of-plane properties
	# This class is the archetype for all other populated material classes.

	def __init__(self, propsDict = None):
		# Collect everything from a dictionary. This first batch are required
		Materials.__init__(self,propsDict)

		# MINIMUM FUNCTIONAL PROPERTY SET
		try:
			self.Name = self.InputDict['name']
			self.E11 = float(self.InputDict['E11'])
			self.E22 = float(self.InputDict['E22'])
			self.E33 = float(self.InputDict['E33'])
			self.Nu12 = float(self.InputDict['Nu12'])
			self.Nu13 = float(self.InputDict['Nu13'])
			self.Nu23 = float(self.InputDict['Nu23'])
			self.G12 = float(self.InputDict['G12'])
			self.G13 = float(self.InputDict['G13'])
			self.G23 = float(self.InputDict['G23'])
			self.Density = float(self.InputDict['Dens'])
			self.CPT = float(self.InputDict['CPT'])
		except KeyError:
			raise KeyError('Check properties input, minimum information not provided')

		# These properties are not required for basic functionality.
		# THERMAL PROPERTIES
		try:
			self.a1 = self.InputDict['a1']
			self.a2 = self.InputDict['a2']
			self.a3 = self.InputDict['a3']
		except KeyError:
			warnings.warn('No CTE included, setting all CTE to zero')
			self.a1 = 0
			self.a2 = 0
			self.a3 = 0

		# Dynamics to be implemented

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

class Plate(Materials):
	# Defines a plate material for use in classic lamination theory.
	# See wiki for where this type is appropriate.

	def __init__(self, propsDict = None):
		# Collect everything from a dictionary. This first batch are required
		Materials.__init__(self, propsDict)

		# MINIMUM FUNCTIONAL PROPERTY SET
		try:
			self.Name = self.InputDict['name']
			self.E11 = float(self.InputDict['E11'])
			self.E22 = float(self.InputDict['E22'])
			self.Nu12 = float(self.InputDict['Nu12'])
			self.G12 = float(self.InputDict['G12'])
			self.Density = float(self.InputDict['Dens'])
			self.CPT = float(self.InputDict['CPT'])
		except KeyError:
			raise KeyError('Check properties input, minimum information not provided')

		# These properties are not required for basic functionality.
		# CTE
		try:
			self.a1 = propsDict['a1']
			self.a2 = propsDict['a2']
		except KeyError:
			warnings.warn('No CTE included, setting all CTE to zero')
			self.a1 = 0
			self.a2 = 0

		# STRESS LIMITS
		try:
			self.f1t = float(self.InputDict['f1t'])
			self.f1c = float(self.InputDict['f1c'])
			self.f2t = float(self.InputDict['f2t'])
			self.f2c = float(self.InputDict['f2c'])
			self.f12s = float(self.InputDict['f12s'])
		except KeyError:
			warnings.warn('No stress limits included, setting all to infinity')
			self.f1t = np.inf
			self.f1c = np.inf
			self.f2t = np.inf
			self.f2c = np.inf
			self.f12s = np.inf

		# STRAIN LIMITS
		try:
			self.e1t = self.InputDict['e1t']
			self.e1c = self.InputDict['e1c']
			self.e2t = self.InputDict['e2t']
			self.e2c = self.InputDict['e2c']
			self.e12s = self.InputDict['e12s']
		except KeyError:
			warnings.warn('No strain limits included, setting all to infinity')
			self.e1t = np.inf
			self.e1c = np.inf
			self.e2t = np.inf
			self.e2c = np.inf
			self.e12s = np.inf

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

class ThickPlate(Plate):
	# Defines a thick plate material for use in first order shear
	# deformation theory. See wiki for details of usage situtaions.

	def __init__(self):
		raise NotImplementedError


class Beam(Materials):

	def __init__(self):
		raise NotImplementedError


class WeakCore(Materials):
	# This is a dummy material that represents a dummy core that adds
	# no strength to the laminate.

	def __init__(self):
		# Weak cores have a name and zero density
		self.Name = 'WeakCore'
		self.Density = 0.0
		self.InputDict = {'name':'WeakCore'}
