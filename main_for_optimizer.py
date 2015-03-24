from laminate_optimizer import *
from composite_materials import *
from laminate_properties import *
import numpy as np

print('LaminateTools v.0.1.0\n (c) Matthew Cowan 2015 \n Released under the Apache v2 License \n')

good = False
while not good:
	print('Enter material properties:')
	E11 = float(raw_input('E11 = '))
	E22 = float(raw_input('E22 = '))
	Nu12 = float(raw_input('nu12 = '))
	G12 = float(raw_input('G12 = '))
	CPT = float(raw_input('CPT = '))

	print('Enter force and moment resultants:')
	Nx = float(raw_input('Nx = '))
	Ny = float(raw_input('Ny = '))
	Nxy = float(raw_input('Nxy = '))
	Mx = float(raw_input('Mx = '))
	My = float(raw_input('My = '))
	Mxy = float(raw_input('Mxy = '))

	print('Enter strain allowables. Unconstrained terms should be set to 1.')
	ex = float(raw_input('ex = '))
	ey = float(raw_input('ey = '))
	exy = float(raw_input('exy = '))

	seed = int(raw_input('Random seed integer = '))

	yn = raw_input('\nInput parameters correct? (y/n)')
	if (yn=='y') or (yn=='Y'):
		good = True



matl = PlateMaterial(name='inputMatl', E11_in=E11, E22_in=E22, Nu12_in=Nu12, G12_in=G12, CPT_in=CPT)
forces = np.matrix([[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]])
allowables = np.matrix([[ex],[ey],[exy],[1],[1],[1]])


result = SimpleOptimizer(material=matl, randSeed=seed, forces=forces, strain_limits=allowables)

print(result['msg'])
print(result['lam'])
