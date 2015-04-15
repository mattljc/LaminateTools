import numpy as np
import materials
import lamination
from copy import deepcopy

class LaminateProperties(object):
	# Superclass containing the basic constructor and essential contract
	# through which the LaminateProperties subclasses operate.

	def __init__(self, lam=None):
		# Constructor wants a Laminate type as its sole input.
		assert isinstance(lam, lamination.Laminate), 'Input not a laminate'
		self.laminate = lam


		# Make global properties common to all types
		self.TotalThickness = 0
		self.TotalAreaDensity = 0
		for ply in self.laminate.PlyStack:
			self.TotalThickness += ply.Thickness
			self.TotalAreaDensity += ply.Material.Density * ply.Thickness

	def __str__(self):
		# Output basic results of the laminate. Ordered keys
		np.set_printoptions(precision=3, linewidth=256, suppress=False)
		output = "== "+self.__class__.__name__+" Properties ==\n"
		keyList = sorted(self.__dict__.keys())

		# Remove keys that aren't useful to report in result string
		while 'specificNT' in keyList: keyList.remove('specificNT')
		while 'laminate' in keyList: keyList.remove('laminate')

		# Build the result string
		for key in keyList:
			output += (key+" = "+str(self.__dict__[key])+"\n")
		return output

	def __repr__(self):
		# Return all the building block representations
		output = '=='+self.__class__.__name__ +'==\n'
		output += repr(self.laminate)
		return output

	def getLaminateProperties():
		# Function generates the requisite laminate properties and
		# stores them in the object.
		raise NotImplementedError

	def calculatePlyStressStrain():
		# Given some sort of global load, calculate the stress and
		# strain in each ply.
		raise NotImplementedError

	def toXML():
		# Outputs useful properties as an XML subtree to be appended to
		# the original XML document.
		raise NotImplementedError

	def toFEA():
		# Outputs material property files for use in NASTRAN, ABAQUS
		# and ANSYS
		raise NotImplementedError

class Continuum(LaminateProperties):
	# Calculates the laminate properties of a thick laminate, including
	# out of plane properties, and shear thinnning effects. This is
	# primarily for use in brick elements.

	# P is defined as the permutation matrix that group inplanes terms.
	# P and its inverse are used extensively in analysis.
	P = np.matrix([ \
	[1,0,0,0,0,0], \
	[0,1,0,0,0,0], \
	[0,0,0,0,0,1], \
	[0,0,0,1,0,0], \
	[0,0,0,0,1,0], \
	[0,0,1,0,0,0]])
	Pinv = P.I

	def __init__(self, lam=None):
		# Call the super-constructor
		LaminateProperties.__init__(lam)
		# This laminate analysis object will have TotalThickness and TotalArealDensity

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isinstance(ply.Material, materials.Continuum), \
			'Wrong material type'

		self.getLaminateProperties()

	def __buildTks(self):
		# Build T-Matrices for each ply and add them to that Ply object
		# This really shouldn't be accessed by itself. Should be
		# accessed as part of property calculations
		for ply in self.laminate.PlyStack:
			m1 = np.cos(np.radians(ply.Orientation))
			m2 = -1* np.sin(np.radians(ply.Orientation))
			n1 = np.sin(np.radians(ply.Orientation))
			n2 = np.cos(np.radians(ply.Orientation))
			#m3 = 0 n3 = 0 p1 = 0 p2 = 0 p3 = 1

			#Simplified transform matrix for quickness and maybe accuracy
			T = np.matrix([ \
			[m1**2, n1**2, 0, 0 , 0 , 2*n1*m1    ], \
			[m2**2, n2**2, 0, 0 , 0 , 2*n2*m2    ], \
			[0    , 0    , 1, 0 , 0 , 0          ], \
			[0    , 0    , 0, n2, m2, 0          ], \
			[0    , 0    , 0, n1, m1, 0          ], \
			[m2*m1, n2*n1, 0, 0 , 0 , (n1*m2+n2*m1)]])

			ply.T = T
			ply.Tinv = T.I

	def getLaminateProperties(self):
		self.__buildTks() # Makes ply transform matrix

		# Build full H-matrices for each ply and store as part of that
		# Ply object. Slice them later.
		for ply in self.laminate.PlyStack:
			ply.GlobalCompliance = ply.Tinv * ply.Material.Compliance * ply.T
			ply.H = self.P * ply.GlobalCompliance.I * self.Pinv

		# Build the A-matrix quadrants. See wiki for details of the math source.
		able = np.matrix( np.zeros((3,3)) )
		baker = np.matrix( np.zeros((3,3)) )
		charlie = np.matrix( np.zeros((3,3)) )
		dog = np.matrix( np.zeros((3,3)) )

		for ply in self.laminate.PlyStack:
			H_II = np.matrix(ply.H[0:3,0:3])
			H_IS = np.matrix(ply.H[0:3,3:])
			H_SI = np.matrix(ply.H[3:,0:3])
			H_SS_inv = np.matrix(ply.H[3:,3:]).I

			able += np.matrix(H_SS_inv * (ply.Thickness / self.TotalThickness))
			baker += np.matrix(H_SS_inv * H_SI * (ply.Thickness / self.TotalThickness))
			charlie += np.matrix(H_IS * H_SS_inv * (ply.Thickness / self.TotalThickness))
			dog += np.matrix((H_II - H_IS * H_SS_inv * H_SI) * (ply.Thickness / self.TotalThickness))

		# Collect the global compliance terms from the A-matrix quadrants
		self.TotalStiffness = np.matrix( np.zeros((6,6)) )
		self.TotalStiffness[0:3,0:3] = dog + charlie * able.I * baker #A_II
		self.TotalStiffness[0:3,3:] = charlie * able.I #A_IS
		self.TotalStiffness[3:,0:3] = able.I * baker #A_SI
		self.TotalStiffness[3:,3:] = able.I #A_SS

		self.TotalCompliance = self.P * self.TotalStiffness.I * self.Pinv
		self.TotalStiffness = self.P * self.TotalStiffness * self.Pinv

		# Use this line if the size of machine epsilon on your machine is problematic
		#self.TotalCompliance = np.matrix(self.TotalCompliance.round(16))

		# Calculate global engineering constants from the global stiffness matrix
		self.Exx = 1 / self.TotalCompliance[0,0]
		self.Eyy = 1 / self.TotalCompliance[1,1]
		self.Ezz = 1 / self.TotalCompliance[2,2]
		self.Gyz = 1 / (2 * self.TotalCompliance[3,3])
		self.Gxz = 1 / (2 * self.TotalCompliance[4,4])
		self.Gxy = 1 / (2 * self.TotalCompliance[5,5])
		self.Nuxy = -self.TotalCompliance[1,0] / self.TotalCompliance[0,0]
		self.Nuxz = -self.TotalCompliance[0,2] / self.TotalCompliance[0,0]
		self.Nuyz = -self.TotalCompliance[1,2] / self.TotalCompliance[1,1]
		self.Etaxs = self.TotalCompliance[0,5] / self.TotalCompliance[0,0]
		self.Etays = self.TotalCompliance[1,5] / self.TotalCompliance[1,1]
		self.Etazs = self.TotalCompliance[2,5] / self.TotalCompliance[2,2]
		self.Etart = self.TotalCompliance[3,4] / (2 * self.TotalCompliance[4,4])
		# NEED: Implementation for 3d CTE

class ThinPlate(LaminateProperties):
	# Calculates the properties of a thin laminate where the
	# out-of-plane properties can be ignored, using classic laminated
	# plate theory.
	def __init__(self, lam=None):
		# Call the super-constructor
		LaminateProperties.__init__(self,lam)

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isinstance(ply.Material, materials.Plate), 'Wrong material type'

		self.getLaminateProperties()

	def getLaminateProperties(self):
		# This analysis assumes a z-axis with zero at the tool surface

		# Initialize A, B and D matrix
		self.A = np.matrix( np.zeros((3,3)) )
		self.B = np.matrix( np.zeros((3,3)) )
		self.D = np.matrix( np.zeros((3,3)) )
		self.specificNT = np.matrix( np.zeros((3,1)) )

		#Initialize zBase
		zLow = -self.TotalThickness / 2

		# Build global compliances ply by ply, then add to the A, B and D matrices
		for ply in self.laminate.PlyStack:
			# Build ply stiffness and CTE in global coordinates
			# Use invariant method to ensure symetry of the final result.
			u1 = ply.Material.U1
			u2 = ply.Material.U2
			u3 = ply.Material.U3
			u4 = ply.Material.U4
			u5 = ply.Material.U5

			c4 = np.cos(4 * np.radians(ply.Orientation))
			c2 = np.cos(2 * np.radians(ply.Orientation))
			c1 = np.cos(1 * np.radians(ply.Orientation))

			s4 = np.sin(4 * np.radians(ply.Orientation))
			s2 = np.sin(2 * np.radians(ply.Orientation))
			s1 = np.sin(1 * np.radians(ply.Orientation))

			Q11 = u1 + u2*c2 + u3*c4
			Q22 = u1 - u2*c2 + u3*c4
			Q12 = u4 - u3*c4
			Q66 = u5 - u3*c4
			Q16 = u2*s2/2 + u3*s4
			Q26 = u2*s2/2 - u3*s4

			ax = ply.Material.a1*c1**2 + ply.Material.a2*s1**2
			ay = ply.Material.a1*s1**2 + ply.Material.a2*c1**2
			axy = (ply.Material.a2-ply.Material.a1)*c1*s1

			ply.GlobalStiffness = np.matrix([
			[Q11, Q12, Q16], \
			[Q12, Q22, Q26], \
			[Q16, Q26, Q66]])

			ply.GlobalCTE = np.matrix([[ax],[ay],[axy]])

			# Build A,B,D matrices
			zUp = zLow + ply.Thickness
			self.A += ply.GlobalStiffness * (zUp - zLow)
			self.B += ply.GlobalStiffness * (zUp**2 - zLow**2) / 2
			self.D += ply.GlobalStiffness * (zUp**3 - zLow**3) / 3
			self.specificNT += ply.GlobalStiffness * ply.GlobalCTE * (zUp - zLow)

			# Increment Z
			zLow = zUp

		self.ABD = np.matrix( np.zeros((6,6)) )
		self.ABD[0:3,0:3] = self.A
		self.ABD[0:3,3:] = self.B
		self.ABD[3:,0:3] = self.B
		self.ABD[3:,3:] = self.D

		# This section generates the effective laminate properties of
		# the plate. Report a soft warining if an asymmetric laminate is
		# given. This is not strictly valid, but effective properties are
		# extremely useful for informed WAGing. Thus we should at least
		# warn the user about what they're doing.
		if self.laminate.Symmetry is False:
			UserWarning('Laminate is not symmetric. Effective properties may not be valid.')

		# Generate properties
		effectiveCompliance = self.A.I
		self.Exx = 1 / (effectiveCompliance[0,0] * self.TotalThickness)
		self.Eyy = 1 / (effectiveCompliance[1,1] * self.TotalThickness)
		self.Gxy = 1 / (effectiveCompliance[2,2] * self.TotalThickness)
		self.Nuxy = - effectiveCompliance[0,1] / effectiveCompliance[0,0]
		self.Etaxs = effectiveCompliance[0,2] / effectiveCompliance[0,0]
		self.Etays = effectiveCompliance[1,2] / effectiveCompliance[1,1]

		effectiveCTE = effectiveCompliance * self.specificNT
		self.ax = effectiveCTE[0,0]
		self.ay = effectiveCTE[1,0]
		self.axy = effectiveCTE[2,0]


	def calculatePlyStressStrain(self, resultants=None):
		# Take force and moment resultants, calculate the global mid-
		# surface strains and curvatures and then iterate through the
		# laminate calculating strains and stresses for each ply.
		# Ply stress/strain is stored in the relevant Ply object.
		# Global strain/curvature is stored in the ThinPlate object.

		# Check type of resultants and enforce matrix shape.
		assert isinstance(resultants, np.matrix), 'resultants is not a matrix'
		assert resultants.shape == (6,1), 'resultants is not 6x1 matrix'

		self.Resultants = resultants
		self.midStrainsCurves = self.ABD.I * self.Resultants
		midStrains = self.midStrainsCurves[:3,0]
		midCurves = self.midStrainsCurves[3:,0]
		# Go through the ply stack, work out the stresses and strains
		# Figuring out z location in exactly the same way as in
		# getLaminateProperties

		zLow = -self.TotalThickness / 2
		for ply in self.laminate.PlyStack:
			# Build the transform matrix
			m = np.cos(np.radians(ply.Orientation))
			n = np.sin(np.radians(ply.Orientation))
			T = np.matrix([\
			[m**2, n**2, 2*m*n],\
			[n**2, m**2, -2*m*n],\
			[-m*n, m*n, m**2-n**2]])

			zUp = zLow + ply.Thickness
			globalStrain = midStrains + zUp*midCurves
			globalStrain[2,0] = globalStrain[2,0]/2 #Convert gamma to epsilon

			ply.Strain = T * globalStrain
			ply.Stress = ply.Material.Compliance.I * ply.Strain
			zLow=zUp

class ThickPlate(LaminateProperties):
	# Calculates the properties of a plate where shear distortion must
	# be included but shear thinning effects can be ignored, such as in
	# a sandwich pannel. This uses first order shear deformation theory.
	def __init__(self, lam=None):
		# Call the super-constructor
		super(Plate,self).__init__(lam)

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isinstance(ply.Material, matl.Continuum)

		self.getLaminateProperties()

# Test Script
if __name__=='__main__':
	t700 ={'name':'T700/2510-UNI', 'E11':18.2e6, 'E22':1.2e6, \
	'G12':0.6e6, 'Nu12':0.3, 'Dens':0.035, 'CPT':0.006}
	t700mat = materials.Plate(t700)

	ply0_2d = lamination.Ply(matl=t700mat, orient=0, thk=0.006)
	lam2d = lamination.Laminate(plyBook=[ply0_2d, deepcopy(ply0_2d)], n_count=1, symmetry=False)
	plate = ThinPlate(lam=lam2d)
	print(str(plate))

	force_moment = np.matrix([0,0,0,10,0,0]).T
	plate.calculatePlyStressStrain(resultants=force_moment)
	count = 1
	for ply in plate.laminate.PlyStack:
		print(id(ply))
		print('Ply #{ct:d} | s1={s1:.3e} s2={s2:.3e} s3={s3:.3e} |'+
		'e1={e1:.3e} e2={e2:.3e} e3={e3:.3e}').format(ct=count,
		s1=ply.Stress[0,0], s2=ply.Stress[1,0], s3=ply.Stress[2,0],
		e1=ply.Strain[0,0], e2=ply.Strain[1,0], e3=ply.Strain[2,0])
		count+=1
