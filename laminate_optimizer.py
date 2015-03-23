import scipy.optimize as opt
import numpy as np
from composite_materials import *
from ply_stack import *
from laminate_properties import *
import time

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
	a = int(lam_details[0])
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
	acceptable = True
	for ct in [0,1,2,3,4,5]:
		acceptable = acceptable and strains[ct]<=strain_limits[ct]

	if acceptable:
		return plate.TotalThickness
	else:
		return np.inf


material = PlateMaterial(name='hw7matl', E11_in=1.85e7, E22_in=1.8e6, Nu12_in=0.3, G12_in=9.3e5, a1_in=0, a2_in=0, ArealDensity_in=0.058, CPT_in=0.006)

ex=1
ey=0.005
exy=0.005
kx = 1.0
ky = 1.0
kxy = 1.0

pressure = 500.00 #psi
radius = 20.0 #in
torque_force = 300000.0 #lb
Nx= pressure * radius
Ny= pressure * radius / 2
Nxy= torque_force / (2*np.pi*radius)
Mx= 0
My= 0
Mxy= 0

strainLimits = np.matrix([[ex],[ey],[exy],[kx],[ky],[kxy]]) # set unconstrained strains to 1, definitely larger than anything that should be calculated
forces = np.matrix([[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]]) # see hand notes for derrivations
initial = [45, 1, 1] # big, massively overbuilt initial guess.

bounds = ((0,90),(1,100),(1,100))

start = time.time()
optimum = opt.differential_evolution(PlateThicknessLimitedStrain, bounds,
args=(material, forces, strainLimits),
popsize = 10,
maxiter=1000,
tol = 0.0001,
mutation=(0.5, 1),
recombination=0.7,
seed=95203,
disp=True,
polish=True)

#bruteAttack = opt.brute(PlateThicknessLimitedStrain, bounds, \
#args=(material, forces, strainLimits), \
#Ns = 20,
#full_output=True, \
#finish = None)
elapsed = time.time() - start

print(optimum.message)
print('a= {a}deg m= {m} n= {n}'.format(a=int(optimum.x[0]), m=int(optimum.x[1]), n=int(optimum.x[2])))
#print(bruteAttack[0])
#print(bruteAttack[1])
print('Time elapsed = {t:.2f}s'.format(t=elapsed))
