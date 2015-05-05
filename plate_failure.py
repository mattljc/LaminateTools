import numpy as np
import materials

def MaxStress(laminate,kind='any'):
	index = list()
	for ply in laminate.PlyStack:
		if isinstance(ply.Material, materials.Plate):
			s1 = ply.Stress[0,0]
			s2 = ply.Stress[1,0]
			s12 = ply.Stress[2,0]

			if (s1>=0):
				f1 = ply.Material.f1t
			else:
				f1 = ply.Material.f1c

			if (s2>=0):
				f2 = ply.Material.f2t
			else:
				f2 = ply.Material.f2c

			f12 = ply.Material.f12s

			k1 = abs(s1/f1)
			k2 = abs(s2/f2)
			k12 = abs(s12/f12)

			ply.MaxStressFail = [k1, k2, k12]
			if kind == 'long': case = k1
			elif kind == 'trans': case = k2
			elif kind == 'shear': case = k12
			else: case = max([k1, k2, k12])

			index.append(case)
	return np.array(index)

def MaxStrain(laminate,kind='any'):
	index = list()
	for ply in laminate.PlyStack:
		if isinstance(ply.Material, materials.Plate):
			eps1 = ply.Strain[0,0]
			eps2 = ply.Strain[1,0]
			eps12 = ply.Strain[2,0]

			e1 = ply.Material.e1t*(e1>=0) + ply.Material.e1c*(e1<0)
			e2 = ply.Material.e2t*(e2>=0) + ply.Material.e2c*(e2<0)
			e12 = ply.Material.e12s

			ply.MaxStrainFail = [eps1/e1, eps2/e2, eps12/e12]
			index.append(max([eps1/e1, eps2/e2, eps12/e12]))
	return np.array(index)

def TsaiHill(laminate,kind='any'):
	index = list()
	for ply in laminate.PlyStack:
		if isinstance(ply.Material, materials.Plate):
			s1 = ply.Stress[0,0]
			s2 = ply.Stress[1,0]
			s12 = ply.Stress[2,0]
			f1 = ply.Material.f1t*(s1>=0) + ply.Material.f1c*(s1<0)
			f2 = ply.Material.f2t*(s2>=0) + ply.Material.f2c*(s2<0)
			f12 = ply.Material.f12s

			ply.TsaiHillFail = (s1/f1)**2 + (s2/f2)**2 + (s12/f12)**2 \
			- s1*s2/(f1**2)
			index.append(ply.TsaiHillFail)
	return np.array(index)

def Hoffman(laminate,kind='any'):
	index = list()
	for ply in laminate.PlyStack:
		if isinstance(ply.Material, materials.Plate):
			s1 = ply.Stress[0,0]
			s2 = ply.Stress[1,0]
			s12 = ply.Stress[2,0]

			f1t = ply.Material.f1t
			f1c = ply.Material.f1c
			f2t = ply.Material.f2t
			f2c = ply.Material.f2c
			f12s = ply.Material.f12s

			ply.HoffmanFail = -s1**2/(f1t*f1c) + s1*s2/(f1t*f1c) \
			- s2**2/(f2t*f2c) + s1*(1/f1t+1/f1c) \
			+ s2*(1/f2t+1/f2c) + (s12/f12s)**2
			index.append(ply.HoffmanFail)
	return np.array(index)
