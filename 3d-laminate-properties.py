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

		self.TotalThickness =

	def buildTk(self):
		# Build a private list containing each directions matrix and its inverse
		# ASSUMPTIONS:
		# 3 Direction is parallel to Z
		T_list = []
		Tinv_list = []

		for ply in laminate:
			orient = ply.Orientation
			m1 = np.cos(np.radians(orient))
			#m2 = np.cos(np.radians(orient+90))
			m2 = -np.sin(np.radians(orient))
			m3 = 0
			#n1 = np.cos(np.radians(90-orient))
			n1 = np.sin(np.radians(orient))
			n2 = np.cos(np.radians(orient))
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

			T_list.append(T)
			Tinv_list.append(T.I)

		self.T = T_list
		self.Tinv = Tinv_list

	def buildHk(self):
		# Build a private list of full H matrices for later slicing
		Sxyz_list = []
		H_list = []
		ct = 0
		for ply in laminate:
			Sk123 = np.matrix([ \
			[1/ply.E11        , -ply.nu12/ply.E22, -ply.nu13/ply.E33, 0            , 0            , 0            ], \
			[-ply.nu12/ply.E11, 1/ply.E22        , -ply.nu32/ply.E33, 0            , 0            , 0            ], \
			[-ply.nu13/ply.E11, -ply.nu23/ply.E22, 1/ply.E33        , 0            , 0            , 0            ], \
			[0                , 0                , 0                , 1/(2*ply.G23), 0            , 0            ], \
			[0                , 0                , 0                , 0            , 1/(2*ply.G13), 0            ], \
			[0                , 0                , 0                , 0            , 0            , 1/(2*ply.G23)]])

			Skxyz = self.Tinv[ct] * Sk123 * self.T[ct]
			Hk = P * Skxyz.I * Pinv
			Sxyz_list.append(Skxyz)
			H_list.append(Hk)

		self.Compliance = Sxyz_list
		self.Hk = H_list

	def getLaminateProperties(self):
