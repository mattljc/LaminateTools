import numpy as np
from composite_materials import *
from ply_stack import *

# beware the size of machine epsilon.... np.finfo(float).eps

class Laminate3D():
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
		# Constructor wants a Laminate type as its sole input.
		if isinstance(lam, Laminate):
			self.laminate = lam
		else:
			raise TypeError('lam is not type Laminate')

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate
			if not isintance(ply.Material, ContinuumMaterial):
				raise TypeError('ply properties are not type ContinuumMaterial')

		# Passed type checks by this point, sum thickness and areal densities
		self.TotalThickness = 0
		self.TotalArealDensity = 0
		for ply in self.laminate.PlyStack:
			self.TotalThickness += ply.Thickness
			self.TotalArealDensity += ply.Material.ArealDensity

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

	def __buildHks(self):
		# Build full H-matrices for each ply and store as part of that Ply object. Slice them later.
		# This really shouldn't be accessed by itself. Should be accessed as part of property calculations

		for ply in self.laminate.PlyStack:
			s11 = 1/ply.Material.E11
			s22 = 1/ply.Material.E22
			s33 = 1/ply.Material.E33
			s12 = -ply.Material.Nu12/ply.Material.E11
			s13 = -ply.Material.Nu13/ply.Material.E11
			s23 = -ply.Material.Nu23/ply.Material.E22
			s66 = 1/(2*ply.Material.G12)
			s55 = 1/(2*ply.Material.G13)
			s44 = 1/(2*ply.Material.G23)

			plyCompliance = np.matrix([ \
			[s11, s12, s13, 0  , 0  , 0  ], \
			[s12, s22, s23, 0  , 0  , 0  ], \
			[s13, s23, s33, 0  , 0  , 0  ], \
			[0  , 0  , 0  , s44, 0  , 0  ], \
			[0  , 0  , 0  , 0  , s55, 0  ], \
			[0  , 0  , 0  , 0  , 0  , s66]])

			ply.GlobalCompliance = ply.Tinv * plyCompliance * ply.T
			ply.H = self.P * ply.GlobalCompliance.I * self.Pinv

	def getLaminateProperties(self):
		# Must call these functions first in this order.
		self.__buildTks() # Makes ply transform matrix
		self.__buildHks() # Makes ply transformed properties matrix

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

	def __str__(self):
		# This string output gets you where you need to go. Just spits out the properties in an easy to read way.
		output = '==3D LAMINATE PROPERTES== \n'
		output += ('E_xx={Exx:.3e}  E_yy={Eyy:.3e}  E_zz={Ezz:.3e}\n').format(Exx=self.Exx, Eyy=self.Eyy, Ezz=self.Ezz)
		output += ('nu_xy={Nuxy:.3g}  nu_xz={Nuxz:.3g}  nu_yz={Nuyz:.3g}\n').format(Nuxy=self.Nuxy, Nuxz=self.Nuxz, Nuyz=self.Nuyz)
		output += ('G_xy={Gxy:.3e}  G_xz={Gxz:.3e}  G_yz={Gyz:.3e}\n').format(Gxy=self.Gxy, Gxz=self.Gxz, Gyz=self.Gyz)
		output += ('eta_xs={Etaxs:g}  eta_ys={Etays:g}  eta_zs={Etazs:g}\n').format(Etaxs=self.Etaxs, Etays=self.Etays, Etazs=self.Etazs)
		output += ('eta_rt={Etart:g}\n').format(Etart=self.Etart)
		return output

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

class Laminate2D():

	def __init__(self, lam=None):

		# Constructor wants a Laminate type as its sole input.
		if isinstance(lam, Laminate):
			self.laminate = lam
		else:
			raise TypeError('lam is not type Laminate')

		# Must check that material definitions are compatable with 3d properties.
		for ply in self.laminate
			if not isintance(ply.Material, PlateMaterial):
				raise TypeError('ply properties are not type ContinuumMaterial')

if __name__ == '__main__':
	# Make a single ply laminate with 0deg orientation to test 3d properties. Output properties should be the same as
	glassUni = RealCompositeMaterial(name='Glass Uni', E11_in=41e9, E22_in=10.4e9, E33_in=10.4e9, Nu12_in=0.28, Nu13_in=0.28, Nu23_in=0.50, G12_in=4.3e9, G13_in=4.3e9, G23_in=3.5e9, ArealDensity_in=1.97, CPT_in=1)
	aPly = Ply(matl=glassUni, orient=90, thk=1)
	thisLam = Laminate(plyBook=[aPly])
	plate = PlateLaminate(lam=thisLam)
	print(plate.verboseString())
