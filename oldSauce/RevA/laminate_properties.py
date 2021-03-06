import numpy as np
import warnings
from composite_materials import *
from ply_stack import *

# beware the size of machine epsilon.... np.finfo(float).eps

class LaminateProperties(object):
	# Superclass containing all the useful overides and the basic constructor

	def __init__(self, lam=None):
		# Constructor wants a Laminate type as its sole input.
		assert isinstance(lam, Laminate)
		self.laminate = lam

		# Make global properties common to all types
		self.TotalThickness = 0
		self.TotalArealDensity = 0
		for ply in self.laminate.PlyStack:
			self.TotalThickness += ply.Thickness
			self.TotalArealDensity += ply.Material.ArealDensity

	def __str__(self):
		# Output basic results of the laminate. Ordered keys
		np.set_printoptions(precision=3, linewidth=128, suppress=False)
		output = "== "+self.__class__.__name__+" Properties ==\n"
		keyList = sorted(self.__dict__.keys())
		while 'specificNT' in keyList: keyList.remove('specificNT')
		while 'laminate' in keyList: keyList.remove('laminate')
		for key in keyList:
			output += (key+" = "+str(self.__dict__[key])+"\n")
		return output

	def __repr__(self):
		# Return all the building block representations
		output = '=='+self.__class__.__name__ +'==\n'
		output += repr(self.laminate)
		return output

	def verbose(self):
		pass

class Laminate3D(LaminateProperties):
	# Calculates the laminate properties of a thick laminate, including out of plane properties.
	# Define permutation matrix P to group inplane terms and its inverse, used extensively in analysis.
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
		super(Laminate3D,self).__init__(lam)
		# This laminate analysis object will have TotalThickness and TotalArealDensity

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isinstance(ply.Material, ContinuumMaterial)

		# Probably not great practice, but construction of this object automatically calls the necessary functions to generate properties.
		self.getLaminateProperties()

	def __buildTks(self):
		# Build T-Matrices for each ply and add them to that Ply object
		# This really shouldn't be accessed by itself. Should be accessed as part of property calculations

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
		# Must call these functions first in this order.
		self.__buildTks() # Makes ply transform matrix

		# Build full H-matrices for each ply and store as part of that Ply object. Slice them later.
		for ply in self.laminate.PlyStack:
			ply.GlobalCompliance = ply.Tinv * ply.Material.Compliance * ply.T
			ply.H = self.P * ply.GlobalCompliance.I * self.Pinv

		# Build the A-matrix quadrants. See wiki for details of the math source
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

	def verboseString(self):
		# This output is a little more intense, spits out everything. Laminate properties and ply by ply compliance matrices
		np.set_printoptions(precision=3)
		output = str(self)
		output += 'Laminate C_xyz =\n'+str(self.TotalCompliance)+'\n'
		output += '\n\n--PLY PROPERTIES--\n'
		ct=1
		for ply in self.laminate.PlyStack:
			output += ('Ply #'+str(ct)+'\n')
			output += (str(ply)+'\n')
			output += ('Ply S_xyz=\n'+str(ply.GlobalCompliance)+'\n')
			output += ('Ply H=\n'+str(ply.H)+'\n')
			output += '----- \n\n'
			ct+=1
		return output

class Laminate2D(LaminateProperties):
	# Calculates the properties of a thin laminate where the out-of-plane properties can be ignored.
	def __init__(self, lam=None):
		# Call the super-constructor
		super(Laminate2D,self).__init__(lam)

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isinstance(ply.Material, PlateMaterial)

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

	def getEffectiveProperties(self):
		# Test that the laminate is symmetric, otherwise these calculations aren't valid
		#assert self.laminate.Symmetry == True

		# Report a soft warining if an asymmetric laminate is given. This is not strictly valid, but effective properties are extremely useful for informed WAGing. Thus we should at least warn the user about what they're doing.
		if (not self.laminate.Symmetry):
			warnings.warn('Laminate is not symmetric. Effective properties may not be valid.', stacklevel=2)

		# Generate properties
		#print(self.A.I)
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

if __name__ == '__main__':
	# 3D Properties: Test Cases: SI UNITS
	# These use metric values for material properties
	# 1 Ply, No rotation (xx=11) and 90 deg rotation (xx=22)
	glassUni3D = ContinuumMaterial(name='Glass Uni', E11_in=41e9, E22_in=10.4e9, E33_in=10.4e9, Nu12_in=0.28, Nu13_in=0.28, Nu23_in=0.50, G12_in=4.3e9, G13_in=4.3e9, G23_in=3.5e9, ArealDensity_in=1.97, CPT_in=1)
	ply_3D_0 = Ply(matl=glassUni3D, orient=0, thk=0.1524e-3)
	ply_3D_90 = Ply(matl=glassUni3D, orient=90, thk=0.1524e-3)
	lam_3D_0 = Laminate(plyBook=[ply_3D_0], n_count=1, symmetry=False)
	lam_3D_90 = Laminate(plyBook=[ply_3D_90], n_count=1, symmetry=False)

	NoRotation3D = Laminate3D(lam=lam_3D_0)
	Rotated3D = Laminate3D(lam=lam_3D_90)

	NoRotated3DTruth = np.isclose(glassUni3D.E11,NoRotation3D.Exx,rtol=0.01) and \
	np.isclose(glassUni3D.E22,NoRotation3D.Eyy,rtol=0.01) and \
	np.isclose(glassUni3D.E33,NoRotation3D.Ezz,rtol=0.01) and \
	np.isclose(glassUni3D.G12,NoRotation3D.Gxy,rtol=0.01) and \
	np.isclose(glassUni3D.G13,NoRotation3D.Gxz,rtol=0.01) and \
	np.isclose(glassUni3D.G23,NoRotation3D.Gyz,rtol=0.01) and \
	np.isclose(glassUni3D.Nu12,NoRotation3D.Nuxy,rtol=0.01) and \
	np.isclose(glassUni3D.Nu13,NoRotation3D.Nuxz,rtol=0.01) and \
	np.isclose(glassUni3D.Nu23,NoRotation3D.Nuyz,rtol=0.01)

	Rotated3DTruth = np.isclose(glassUni3D.E11,Rotated3D.Eyy,rtol=0.01) and \
	np.isclose(glassUni3D.E22,Rotated3D.Exx,rtol=0.01) and \
	np.isclose(glassUni3D.E33,Rotated3D.Ezz,rtol=0.01) and \
	np.isclose(glassUni3D.G12,Rotated3D.Gxy,rtol=0.01) and \
	np.isclose(glassUni3D.G13,Rotated3D.Gyz,rtol=0.01) and \
	np.isclose(glassUni3D.G23,Rotated3D.Gxz,rtol=0.01) and \
	np.isclose((glassUni3D.Nu12*glassUni3D.E22/glassUni3D.E11),Rotated3D.Nuxy,rtol=0.01) and \
	np.isclose(glassUni3D.Nu13,Rotated3D.Nuyz,rtol=0.01) and \
	np.isclose(glassUni3D.Nu23,Rotated3D.Nuxz,rtol=0.01)
	# Note that Nuxy is rotated so that Nuxy = Nu21 != Nu12, thus the comparison is to the Nu21 value calculated based on the symmetry of the compliance matrix. See wiki.

	print('3D LAMINATE TEST RESULTS\n+++++++++++++++++++++++++++++++++++++++')
	#print(NoRotation3D)
	#print(Rotated3D.verboseString())
	#print(Rotated3D)
	print('3D 1-ply, 0deg Pass? '+str(NoRotated3DTruth))
	print('3D 1-ply, 90deg Pass? '+str(Rotated3DTruth))

	############################################################################
	# 2D Properties: Test Cases: SI UNITS
	# 1 Ply, No rotation (xx=11) and 90 deg rotation (xx=22)
	glassUni2D = PlateMaterial(name='Glass Uni Plate', E11_in=41e9, E22_in=10.4e9, Nu12_in=0.28, G12_in=4.3e9, a1_in= 7e-6, a2_in= 26e-6,ArealDensity_in=1.97, CPT_in=1)
	ply_2D_0 = Ply(matl=glassUni2D, orient=0, thk=0.1524e-3)
	ply_2D_90 = Ply(matl=glassUni2D, orient=90, thk=0.1524e-3)
	lam_2D_0_S = Laminate(plyBook=[ply_2D_0], n_count=1, symmetry=True)
	lam_2D_90_S = Laminate(plyBook=[ply_2D_90], n_count=1, symmetry=True)
	# Need to generate a symetric 3D composite for comparison
	lam_3D_0_S = Laminate(plyBook=[ply_3D_0], n_count=1, symmetry=True)
	lam_3D_90_S = Laminate(plyBook=[ply_3D_90], n_count=1, symmetry=True)

	NoRotation2D = Laminate2D(lam=lam_2D_0_S)
	NoRotation2D3D = Laminate3D(lam=lam_3D_0_S)
	NoRotation2D.getEffectiveProperties()

	Rotated2D = Laminate2D(lam=lam_2D_90_S)
	Rotated2D3D = Laminate3D(lam=lam_3D_90_S)
	Rotated2D.getEffectiveProperties()

	NoRotated2DTruth = np.isclose(NoRotation2D.Exx,NoRotation2D3D.Exx,rtol=0.01)
	NoRotated2DTruth = NoRotated2DTruth and np.isclose(NoRotation2D.Eyy,NoRotation2D3D.Eyy,rtol=0.01)
	NoRotated2DTruth = NoRotated2DTruth and np.isclose(NoRotation2D.Gxy,NoRotation2D3D.Gxy,rtol=0.01)
	NoRotated2DTruth = NoRotated2DTruth and np.isclose(NoRotation2D.Nuxy,NoRotation2D3D.Nuxy,rtol=0.01)

	Rotated2DTruth = np.isclose(Rotated2D.Exx,Rotated2D3D.Exx,rtol=0.01)
	Rotated2DTruth = Rotated2DTruth and np.isclose(Rotated2D.Eyy,Rotated2D3D.Eyy,rtol=0.01)
	Rotated2DTruth = Rotated2DTruth and np.isclose(Rotated2D.Gxy,Rotated2D3D.Gxy,rtol=0.01)
	Rotated2DTruth = Rotated2DTruth and np.isclose(Rotated2D.Nuxy,Rotated2D3D.Nuxy,rtol=0.01)

	print('\n2D LAMINATE TEST RESULTS\n+++++++++++++++++++++++++++++++++++++++')
	#print(NoRotation2D)
	#print(Rotated2D)
	print('2D 1-ply, 0deg Pass? '+str(NoRotated2DTruth))
	print('2D 1-ply, 90deg Pass? '+str(Rotated2DTruth))
