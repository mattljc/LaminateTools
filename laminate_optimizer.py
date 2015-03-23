import scipy.optimize as opt
import numpy as np
from composite_materials import *
from ply_stack import *
from laminate_properties import *

def PlateThicknessLimitedStrain(lam_details, material, forces, strain_limits):
	# Fitness returns the thickness of a laminate having analyzed its properties and calculated overall strains
	# If the calculated strains exceed strain_limits the thickness is returned as infinite
	# Assumes a laminate of the form [a_m/90_n/0_n/90_n/-a_m]
	# INPUT DEFINITIONS:
	# lam_details: the laminate definition as a list or tupple in the following form [a, n, m]
	# material: a PlateMaterial object
	# forces: a list of thickness integrated forces and moments in the form [Nx, Ny, Nxy, Mx, My, Mxy]
	# strain_limits: a list of global strain limits in the form [ex, ey, exy, kx, ky, kxy]

	# Cast inputs as appropriate types
	a = lam_details[0]
	m = int(lam_details[1])
	n = int(lam_details[2])
	#print('[{a:.2g}_{m} / 0_{n} / 90_{n} / 0_{n} / {negA:.2g}_{m} ]'.format(a=a, m=m, n=n, negA=a*-1))

	# Build the plies
	plyA = Ply(matl=material, orient=a)
	plyNegA = Ply(matl=material, orient=-a)
	ply0 = Ply(matl=material, orient=0)
	ply90 = Ply(matl=material, orient=90)

	# Build the plybook, laminate and analyses
	plybook = [plyA]*m + [ply0]*n + [ply90]*n + [ply0]*n + [plyNegA]*m
	laminate = Laminate(plyBook=plybook)
	plate = Laminate2D(lam=laminate)

	# Take output a,b,d matrices from analysis and build an augmented ABD matrix to determine strains
	augmentedABD = np.matrix( np.zeros((6,6)) )
	augmentedABD[0:3,0:3] = plate.A
	augmentedABD[0:3,3:] = plate.B
	augmentedABD[3:,0:3] = plate.B
	augmentedABD[3:,3:] = plate.D

	# Calculate global strains and compare to limits
	strains = augmentedABD.I * forces
	#print(strains.T)
	acceptable = np.allclose(strains,strain_limits)
	if acceptable:
		return plate.TotalThickness
	else:
		return 1e9


material = PlateMaterial(name='hw7matl', E11_in=20e6, E22_in=1e6, Nu12_in=0.3, G12_in=0.7e6, a1_in=0, a2_in=0, ArealDensity_in=0.058, CPT_in=0.006)

ex=2e-5
ey=2e-5
exy=2e-5
kx = 1.0
ky = 1.0
kxy = 1.0

pressure = 200.00 #psi
radius = 20.0 #in
torque_force = 300 #lb
Nx= pressure * radius
Ny= pressure * radius / 2
Nxy= torque_force / (2*np.pi)
Mx= 0
My= 0
Mxy= 0

strainLimits = np.matrix([[ex],[ey],[exy],[kx],[ky],[kxy]]) # set unconstrained strains to 1, definitely larger than anything that should be calculated
forces = np.matrix([[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]]) # see hand notes for derrivations
initial = [45, 1, 1] # big, massively overbuilt initial guess.

bounds = ((-90,90),(1,20),(1,20))
#optimum = opt.minimize(PlateThicknessLimitedStrain, initial, args=(material, forces, strainLimits), method='Nelder-Mead', options={'disp': True,'ftol':1e-8})
optimum = opt.differential_evolution(PlateThicknessLimitedStrain,bounds,args=(material, forces, strainLimits),maxiter=10)
#bruteAttack =

print(optimum.message)
print(optimum.x)
