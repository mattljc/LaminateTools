import numpy as np
from composite_materials import *
from ply_stack import *

# beware the size of machine epsilon.... np.finfo(float).eps

class Laminate():
	# Superclass containing all the useful overides and the basic constructor

	def __init__(self):
		# Constructor wants a Laminate type as its sole input.
		assert isinstance(lam, Laminate)

		# Make global properties common to all types
		self.TotalThickness = 0
		self.TotalArealDensity = 0
		for ply in self.laminate.PlyStack:
			self.TotalThickness += ply.Thickness
			self.TotalArealDensity += ply.Material.ArealDensity

	def __str__(self):
		# Output basic results of the laminate. Ordered keys
		output = "== "+self.__class__.__name__+" Properties =="
		for key in sorted(self.__dict__.keys())
			output += ("{property} = {val:.3e}\n").format(property=key, val=self.__dict__[key])
		return output

	def __repr__(self):
		# Return all the building block representations
		output = self.__class__.__name__ + '\n'
		output += repr(self.Laminate)
		return output

	def verbose(self):


class Laminate3D(Laminate):
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
		self = Laminate(lam)
		# This laminate analysis object will have TotalThickness and TotalArealDensity

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate:
			assert isintance(ply.Material, ContinuumMaterial)

		# Probably not great practice, but construction of this object automatically calls the necessary functions to generate properties.
		self.getLaminateProperties()

	def __buildTks(self):
		# Build T-Matrices for each ply and add them to that Ply object
		# This really shouldn't be accessed by itself. Should be accessed as part of property calculations

		for ply in self.laminate.PlyStack:
			m1 = np.cos(np.radians(ply.Orientation))
			#m2 = np.cos(np.radians(orient+90))
			m2 = -np.sin(np.radians(ply.Orientation))
			m3 = 0
			#n1 = np.cos(np.radians(90-orient))
			n1 = np.sin(np.radians(ply.Orientation))
			n2 = np.cos(np.radians(ply.Orientation))
			n3 = 0
			p1 = 0
			p2 = 0
			p3 = 1

			#Full Transform matrix for reference
			#T = np.matrix([[m1**2, n1**2, p1**2, 2*n1*p1, 2*m1*p1, 2*n1*m1], \
			#               [m2**2, n2**2, p2**2, 2*n2*p2, 2*m2*p2, 2*n2*m2], \
			#               [m3**2, n3**2, p3**2, 2*n3*p3, 2*m3*p3, 2*n3*m3], \
			#               [m2*m3, n2*n3, p2*p3, n2*p3+n3*p2, m3*p2+m2*p3, n2*m3+n3*m2], \
			#               [m1*m3, n1*n3, p1*p3, n1*p3+n3*p1, m1*p3+m3*p1, n1*m3+n3*m1], \
			#               [m2*m1, n2*n1, p2*p1, n1*p2+n2*p1, m1*p2+m2*p1, n1*m2+n2*m1]])

			#Simplified transform matrix for quickness and maybe accuracy
			T = np.matrix([ \
			[m1**2, n1**2, 0, 0 , 0 , 2*n1*m1    ], \
			[m2**2, n2**2, 0, 0 , 0 , 2*n2*m2    ], \
			[0    , 0    , 1, 0 , 0 , 0          ], \
			[0    , 0    , 0, n2, m2, 0          ], \
			[0    , 0    , 0, n1, m1, 0          ], \
			[m2*m1, n2*n1, 0, 0 , 0 , n1*m2+n2*m1]])

			ply.T = T
			ply.Tinv = T.I

	def getLaminateProperties(self):
		# Must call these functions first in this order.
		self.__buildTks() # Makes ply transform matrix

		# Build full H-matrices for each ply and store as part of that Ply object. Slice them later.
		for ply in self.laminate.PlyStack:
			ply.GlobalCompliance = ply.Tinv * ply.Material.Compliance * ply.T
			ply.H = self.P * ply.GlobalCompliance.I * self.Pinv

		# Build the A-matrix quadrants. See wiki for details of the math
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
		self.TotalCompliance = np.matrix( np.zeros((6,6)) )
		self.TotalCompliance[0:3,0:3] = dog + charlie * able.I * baker #A_II
		self.TotalCompliance[0:3,3:] = charlie * able.I #A_IS
		self.TotalCompliance[3:,0:3] = able.I * baker #A_SI
		self.TotalCompliance[3:,3:] = able.I #A_SS

		self.TotalCompliance = self.P * self.TotalCompliance.I * self.Pinv
		# Use this line if the size of machine epsilon on your machine is problematic
		#self.TotalCompliance = np.matrix(self.TotalCompliance.round(16))

		# Calculate global engineering constants from the global stiffness matrix
		self.Exx = 1 / self.TotalCompliance[0,0]
		self.Eyy = 1 / self.TotalCompliance[1,1]
		self.Ezz = 1 / self.TotalCompliance[2,2]
		self.Gyz = 1 / (2 * self.TotalCompliance[3,3])
		self.Gxz = 1 / (2 * self.TotalCompliance[4,4])
		self.Gxy = 1 / (2 * self.TotalCompliance[5,5])
		self.Nuxy = self.TotalCompliance[0,1] * (- self.Exx)
		self.Nuxz = self.TotalCompliance[0,2] * (- self.Exx)
		self.Nuyz = self.TotalCompliance[1,2] * (- self.Eyy)
		self.Etaxs = self.TotalCompliance[0,5] * self.Exx
		self.Etays = self.TotalCompliance[1,5] * self.Eyy
		self.Etazs = self.TotalCompliance[2,5] * self.Ezz
		self.Etart = self.TotalCompliance[3,4] * 2 * self.Gxz

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

class Laminate2D(Laminate):

	def __init__(self, lam=None):
		# Call the super-constructor
		self = Laminate(lam)

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate.PlyStack:
			assert isintance(ply.Material, PlateMaterial)

	def getLaminateProperties(self):
		# This analysis assumes a z-axis with zero at the tool surface
		# Figure out the midplane ordinate
		midplane = self.TotalThickness / 2

		# Initialize A, B and D matrix
		self.A = np.matrix( np.zeros((3,3)) )
		self.B = np.matrix( np.zeros((3,3)) )
		self.D = np.matrix( np.zeros((3,3)) )

		#Initialize zBase
		zLow = midplane

		# Build global compliances ply by ply, then add to the A, B and D matrices
		for ply in self.Laminate.PlyStack

			zUp = abs(zLow - ply.Thickness)

			c2 = np.cos(np.radians(2* ply.Orientation))
			c4 = np.cos(np.radians(4* ply.Orientation))
			s2 = np.sin(np.radians(2* ply.Orientation))
			s4 = np.sin(np.radians(4* ply.Orientation))

			q11 = ply.Material.U1 + c2*ply.Material.U2 + c4*ply.Material.U3
			q12 = ply.Material.U4 + c4*ply.Material.U3
			q16 = s2*ply.Material.U2/2 + s4*ply.Material.U3
			q22 = ply.Material.U1 - c2*ply.Material.U2 + c4*ply.Material.U3
			q26 = s2*ply.Material.U2/2 - s4*ply.Material.U3
			q66 =

			ply.GlobalCompliance = np.Matrix([
			[q11, q12, q16], \
			[q12, q22, q26], \
			[q16, q26, q66]])

			self.A += ply.GlobalCompliance * (zUp - zLow)
			self.B += ply.GlobalCompliance * (zUp**2 - zLow**2) / 2
			self.D += ply.GlobalCompliance * (zUp**3 - zLow**3) / 3

			zLow = zUp

	def getEffectiveProperties(self):
		# Test that the laminate is symmetric, otherwise these calculations aren't valid
		assert self.laminate.Symmetry == True

		# Generate properties


if __name__ == '__main__':
	# Make a single ply laminate with 0deg orientation to test 3d properties. Output properties should be the same as
	glassUni = RealCompositeMaterial(name='Glass Uni', E11_in=41e9, E22_in=10.4e9, E33_in=10.4e9, Nu12_in=0.28, Nu13_in=0.28, Nu23_in=0.50, G12_in=4.3e9, G13_in=4.3e9, G23_in=3.5e9, ArealDensity_in=1.97, CPT_in=1)
	aPly = Ply(matl=glassUni, orient=90, thk=1)
	thisLam = Laminate(plyBook=[aPly])
	thisPlate = Laminate3D(lam=thisLam)
	print(thisPlate.verboseString())
