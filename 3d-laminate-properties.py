import numpy as np
from ply_stack import *

#beware the size of machine epsilon.... np.finfo(float).eps

class PlateLaminate():

	#Define permutation matrix P to group inplane terms and its inverse
	P = np.matrix([ \
	[1,0,0,0,0,0,0], \
	[0,1,0,0,0,0,0], \
	[0,0,0,0,0,0,1], \
	[0,0,0,1,0,0,0], \
	[0,0,0,0,1,0,0], \
	[0,0,1,0,0,0,0]])
	Pinv = P.I

	def __init__(self, lam=None):
		if isinstance(lam, Laminate):
			self.__laminate = lam
		else:
			raise TypeError('lam is not type Laminate')

		self.TotalThickness = 0
		for ply in self.laminate:
			self.TotalThickness += ply.Thickness

	def buildTks(self):
		# Build T-Matrices for each ply and add them to that Ply object
		# ASSUMPTIONS:
		# 3 Direction is parallel to Z

		for ply in self.Laminate:
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
			[0    , 0    , 1, n2, m2, 0          ], \
			[0    , 0    , 0, n1, m1, 0          ], \
			[m2*m1, n2*n1, 0, 0 , 0 , n1*m2+n2*m1]])

			ply.T = T
			ply.Tinv = T.I

	def buildHks(self):
		# Build full H-matrices for each ply and store as part of that Ply object. Slice them later.
		Sxyz_list = []
		H_list = []
		ct = 0
		for ply in self.Laminate:
			Sk123 = np.matrix([ \
			[1/ply.E11        , -ply.nu12/ply.E22, -ply.nu13/ply.E33, 0            , 0            , 0            ], \
			[-ply.nu12/ply.E11, 1/ply.E22        , -ply.nu32/ply.E33, 0            , 0            , 0            ], \
			[-ply.nu13/ply.E11, -ply.nu23/ply.E22, 1/ply.E33        , 0            , 0            , 0            ], \
			[0                , 0                , 0                , 1/(2*ply.G23), 0            , 0            ], \
			[0                , 0                , 0                , 0            , 1/(2*ply.G13), 0            ], \
			[0                , 0                , 0                , 0            , 0            , 1/(2*ply.G23)]])

			ply.Compliance = ply.Tinv * Sk123 * ply.T
			ply.H = P * ply.Compliance.I * Pinv

	def getLaminateProperties(self):
		
