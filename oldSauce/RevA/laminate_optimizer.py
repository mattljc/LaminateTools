import scipy.optimize as opt
import numpy as np
from composite_materials import *
from ply_stack import *
from laminate_properties import *
import time

def PlateThicknessLimitedStrain(lam_details, material, forces, strain_limits):
	# Fitness returns the thickness of a laminate having analyzed its properties and calculated overall strains
	# If the calculated strains exceed strain_limits the thickness is returned as infinite
	# Assumes a laminate of the form [a_m/-a_m/0_n/__90_n]s
	# INPUT DEFINITIONS:
	# lam_details: the laminate definition as a list or tupple in the following form [a, n, m]
	# material: a PlateMaterial object
	# forces: a list of thickness integrated forces and moments in the form [Nx, Ny, Nxy, Mx, My, Mxy]
	# strain_limits: a list of global strain limits in the form [ex, ey, exy, kx, ky, kxy]

	# Cast inputs as appropriate types
	a = lam_details[0]
	m = int(lam_details[1])
	n = int(lam_details[2])

	laminate = LaminateBuilder(a,m,n,material)
	plate = Laminate2D(lam=laminate)

	# Calculate global strains and compare to limits
	strains = plate.ABD.I * forces
	#print(strains.T)
	acceptable = True
	for ct in [0,1,2,3,4,5]:
		acceptable = acceptable and strains[ct]<=strain_limits[ct]

	if acceptable:
		thk=plate.TotalThickness
	else:
		thk=np.inf

	#print(('[a={ao} m={mo} n={no}] -- [ey={eyo:.3e} exy={exyo:.3e}] -- thk={to}').format(ao=a,mo=m,no=n,eyo=strains[1,0],exyo=strains[2,0],to=thk))
	return thk

def LaminateBuilder(a=45,m=10,n=10, material=None):
	# Build the plies
	plyA = Ply(matl=material, orient=a)
	plyNegA = Ply(matl=material, orient=-a)
	ply0 = Ply(matl=material, orient=0)
	ply90 = Ply(matl=material, orient=90)

	# Build the plybook, laminate and analyses
	plybook = [plyA]*m + [plyNegA]*m + [ply0]*n + [ply90]*n + [ply0]*n + [plyNegA]*m + [plyA]*m
	laminate = Laminate(plyBook=plybook)

	return laminate

def SimpleOptimizer(material=None, bounds=((0,90),(1,100),(1,100)), randSeed=599, forces=np.matrix([[1],[1],[1],[1],[1],[1]]), strain_limits=np.matrix([[1],[1],[1],[1],[1],[1]])):
	# Simple optimizer is the entry point for this optimization system.
	optimum = opt.differential_evolution(PlateThicknessLimitedStrain, bounds,
	args=(material, forces, strain_limits),
	popsize = 10,
	maxiter=5000,
	tol = 0.00001,
	mutation=(0.5, 1),
	seed = randSeed,
	recombination=0.7,
	disp=True,
	polish=False)

	finalLam = LaminateBuilder(a=int(optimum.x[0]), m=int(optimum.x[1]), n=int(optimum.x[2]), material=material)
	finalPlate = Laminate2D(lam=finalLam)
	finalPlate.getEffectiveProperties()
	finalStrains = finalPlate.ABD.I*forces

	print('\n'+optimum.message)

	output = 'a= {a:.3f}deg m= {m} n= {n}'.format(a=optimum.x[0], m=int(optimum.x[1]), n=int(optimum.x[2]))
	output += '\nNx={nx:.3e} Ny={ny:.3e} Nxy={nxy:.3e}'.format(nx=forces[0,0], ny=forces[1,0], nxy=forces[2,0])
	output += '\nex={ex:.3e} ey={ey:.3e} exy={exy:.3e}'.format(ex=finalStrains[0,0], ey=finalStrains[1,0], exy=finalStrains[2,0])
	output += '\nOptimized Laminate Properties'
	output += '\n'+str(finalPlate)

	return {'msg':output,'lam':finalLam}

if __name__ == '__main__':
	material = PlateMaterial(name='hw7matl', E11_in=1.85e7, E22_in=1.8e6, Nu12_in=0.3, G12_in=9.3e5, a1_in=0, a2_in=0, ArealDensity_in=0.058, CPT_in=0.006)

	ex=1
	ey=0.005
	exy=0.005
	kx = 1.0
	ky = 1.0
	kxy = 1.0

	pressure = 500.00 #psi
	radius = 10.0 #in
	torque_force = 3e5 #lb
	Nx= pressure * radius
	Ny= pressure * radius / 2
	Nxy= torque_force / (2*np.pi*radius)
	Mx= 0
	My= 0
	Mxy= 0

	strainLimits = np.matrix([[ex],[ey],[exy],[kx],[ky],[kxy]]) # set unconstrained strains to 1, definitely larger than anything that should be calculated
	forces = np.matrix([[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]]) # see hand notes for derrivations
	initial = [45, 100, 100] # big, massively overbuilt initial guess.

	bounds = ((0,90),(1,100),(1,100))

	optimum = opt.differential_evolution(PlateThicknessLimitedStrain, bounds,
	args=(material, forces, strainLimits),
	popsize = 10,
	maxiter=5000,
	tol = 0.00001,
	mutation=(0.5, 1),
	seed = 692962,
	recombination=0.7,
	disp=True,
	polish=False)


	finalLam = LaminateBuilder(a=int(optimum.x[0]), m=int(optimum.x[1]), n=int(optimum.x[2]), material=material)
	finalPlate = Laminate2D(lam=finalLam)
	finalPlate.getEffectiveProperties()
	finalStrains = finalPlate.ABD.I*forces

	print('\n'+optimum.message)
	print('a= {a:.3f}deg m= {m} n= {n}'.format(a=optimum.x[0], m=int(optimum.x[1]), n=int(optimum.x[2])))
	print('Nx={nx:.3e} Ny={ny:.3e} Nxy={nxy:.3e}'.format(nx=Nx, ny=Ny, nxy=Nxy))
	print('ex={ex:.3e} ey={ey:.3e} exy={exy:.3e}'.format(ex=finalStrains[0,0], ey=finalStrains[1,0], exy=finalStrains[2,0]))
	print('\nOptimized Laminate Properties')
	print(finalPlate)
