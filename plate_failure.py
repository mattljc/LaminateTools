def MaxStress(laminate):
	for ply in laminate.PlyBook:
		s1 = ply.Stress[0]
		s2 = ply.Stress[1]
		s12 = ply.Stress[2]

		f1t = ply.Material.f1t
		f1c = ply.Material.f1c
		f2t = ply.Material.f2t
		f2c = ply.Material.f2c
		f12s = ply.Material.f12s

		index = [s1/f1t, s1/f1c, s2/f2t, s2/f2c, s12/f12s]
		ply.MaxStressFail = max(index)

def MaxStrain(laminate):
	for ply in laminate.PlyBook:
		e1 = ply.Strain[0]
		e2 = ply.Strain[1]
		e12 = ply.Strain[2]

		e1t = ply.Material.f1t
		e1c = ply.Material.f1c
		e2t = ply.Material.f2t
		e2c = ply.Material.f2c
		e12s = ply.Material.f12s

		index = [e1/e1t, e1/e1c, e2/e2t, e2/e2c, e12/e12s]
		ply.MaxStrainFail = max(index)

def TsaiHill(laminate):
	for ply in laminate.PlyBook:
		s1 = ply.Stress[0]
		s2 = ply.Stress[1]
		s12 = ply.Stress[2]
		f1 = ply.Material.f1t*(s1>=0) + ply.Material.f1c*(s1<0)
		f2 = ply.Material.f2t*(s2>=0) + ply.Material.f2c*(s2<0)
		f12 = ply.Material.f12s

		ply.TsaiHillFail = (s1/f1)**2 + (s2/f2)**2 + (s12/f12)**2 \
		- s1*s2/(f1**2)

def Hoffman(laminate):
	for ply in laminate.PlyBook:
		s1 = ply.Stress[0]
		s2 = ply.Stress[1]
		s12 = ply.Stress[2]

		f1t = ply.Material.f1t
		f1c = ply.Material.f1c
		f2t = ply.Material.f2t
		f2c = ply.Material.f2c
		f12s = ply.Material.f12s

		ply.HoffmanFail = -s1**2/(f1t*f1c) + s1*s2/(f1t*f1c) \
		- s2**2/(f2t*f2c) + s1*(f1c+f1t)/(f1c*f1t) \
		+ s2*(f2c+f2t)/(f2c*f2t) + (s12/f12s)**2
