import numpy as np
from composite_materials import *
from ply_stack import *
from laminate_properties import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt


def LaminateBuilder(a=45,m=10,n=10,material=None):
	# Build the plies
	plyA = Ply(matl=material, orient=a)
	plyNegA = Ply(matl=material, orient=-a)
	ply0 = Ply(matl=material, orient=0)
	ply90 = Ply(matl=material, orient=90)

	# Build the plybook, laminate and analyses
	plybook = plybook = [plyA]*m + [plyNegA]*m + [ply0]*n + [ply90]*n + [ply0]*n + [plyNegA]*m + [plyA]*m
	laminate = Laminate(plyBook=plybook)

	return laminate

# Make the ranges and mesh grids
# Avoid the curse of the off-by-one error
aRange = np.arange(0,91,5)
mRange = np.arange(1,52,2)
nRange = np.arange(1,52,2)
aBaseline = 46.153
mBaseline = 8
nBaseline = 2
a_amMesh, m_amMesh = np.meshgrid(aRange, mRange)
a_anMesh, n_anMesh = np.meshgrid(aRange, nRange)
m_mnMesh, n_mnMesh = np.meshgrid(mRange, nRange)

# Initialize collector arrays for strains
ey_am = np.zeros(a_amMesh.shape)
exy_am = np.zeros(a_amMesh.shape)
ey_an = np.zeros(a_anMesh.shape)
exy_an = np.zeros(a_anMesh.shape)
ey_mn = np.zeros(m_mnMesh.shape)
exy_mn = np.zeros(m_mnMesh.shape)

# Calculate integrated forces and moments. Make material
pressure = 500.00 #psi
radius = 10.0 #in
torque_force = 300000.0 #lb
Nx= pressure * radius
Ny= pressure * radius / 2
Nxy= torque_force / (2*np.pi*radius)
Mx= 0
My= 0
Mxy= 0
material = PlateMaterial(name='hw7matl', E11_in=1.85e7, E22_in=1.8e6, Nu12_in=0.3, G12_in=9.3e5, a1_in=0, a2_in=0, ArealDensity_in=0.058, CPT_in=0.006)
forces = np.matrix([[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]])

# Start some monstrous iteration
# a-m plane
for row in np.arange(0,a_amMesh.shape[0],1): #Rows
	for col in np.arange(0,a_amMesh.shape[1],1): #Columns
		lam = LaminateBuilder(a=a_amMesh[row,col], m=m_amMesh[row,col], n=nBaseline, material=material)
		props = Laminate2D(lam=lam)
		strains = props.ABD.I * forces
		ey_am[row,col] = strains[1,0]
		exy_am[row,col] = strains[2,0]

#a-n plane
for row in np.arange(0,a_anMesh.shape[0],1): #Rows
	for col in np.arange(0,a_anMesh.shape[1],1): #Columns
		lam = LaminateBuilder(a=a_anMesh[row,col], m=mBaseline, n=n_anMesh[row,col], material=material)
		props = Laminate2D(lam=lam)
		strains = props.ABD.I * forces
		ey_an[row,col] = strains[1,0]
		exy_an[row,col] = strains[2,0]

#m-n plane
for row in np.arange(0,m_mnMesh.shape[0],1): #Rows
	for col in np.arange(0,m_mnMesh.shape[1],1): #Columns
		lam = LaminateBuilder(a=aBaseline, m=m_mnMesh[row,col], n=n_mnMesh[row,col], material=material)
		props = Laminate2D(lam=lam)
		strains = props.ABD.I * forces
		ey_mn[row,col] = strains[1,0]
		exy_mn[row,col] = strains[2,0]

#print(ey_am)
#print(ey_an)
#print(ey_mn)

lvs = [.005, .004, .003, .002, .001]

# Plot a-m plane
plt.figure(1)
am_ey_contour = plt.contour(a_amMesh,m_amMesh,ey_am, levels=lvs)
plt.plot(aBaseline,mBaseline,'ko')
plt.clabel(am_ey_contour, inline=1, fontsize=10)
#plt.title('a-m plane\nLongitudinal Strains, Constant n={const}'.format(const=mBaseline))
plt.axis([0,90,0,50])
plt.xlabel('a [deg]')
plt.ylabel('m [# plies]')
pylab.savefig('amTrans.pdf', bbox_inches='tight')
############
plt.figure(2)
am_exy_contour = plt.contour(a_amMesh,m_amMesh,exy_am, levels=lvs)
plt.plot(aBaseline,mBaseline,'ko')
plt.clabel(am_exy_contour, inline=1, fontsize=10)
#plt.title('a-m plane\nShear Strains, Constant n={const}'.format(const=mBaseline))
plt.axis([0,90,0,50])
plt.xlabel('a [deg]')
plt.ylabel('m [# plies]')
pylab.savefig('amShear.pdf', bbox_inches='tight')

# Plot a-n plane
plt.figure(3)
an_ey_contour = plt.contour(a_anMesh,n_anMesh,ey_an, levels=lvs)
plt.plot(aBaseline,nBaseline,'ko')
plt.clabel(an_ey_contour, inline=1, fontsize=10)
#plt.title('a-n plane\nLongitudinal Strains, Constant m={const}'.format(const=nBaseline))
plt.axis([0,90,0,50])
plt.xlabel('a [deg]')
plt.ylabel('n [# plies]')
pylab.savefig('anTrans.pdf', bbox_inches='tight')
#############
plt.figure(4)
an_exy_contour = plt.contour(a_anMesh,n_anMesh,exy_an, levels=lvs)
plt.plot(aBaseline,nBaseline,'ko')
plt.clabel(an_exy_contour, inline=1, fontsize=10)
#plt.title('a-n plane\nShear Strains, Constant m={const}'.format(const=nBaseline))
plt.axis([0,90,0,50])
plt.xlabel('a [deg]')
plt.ylabel('n [# plies]')
pylab.savefig('anShear.pdf', bbox_inches='tight')

# Plot m-n plane
plt.figure(5)
mn_ey_contour = plt.contour(m_mnMesh,n_mnMesh,ey_mn, levels=lvs)
plt.plot(mBaseline,nBaseline,'ko')
plt.clabel(mn_ey_contour, inline=1, fontsize=10)
#plt.title('m-n plane\nLongitudinal Strains, Constant a={const}'.format(const=aBaseline))
plt.axis([0,50,0,50])
plt.xlabel('m [# plies]')
plt.ylabel('n [# plies]')
pylab.savefig('mnTrans.pdf', bbox_inches='tight')
#############
plt.figure(6)
mn_exy_contour = plt.contour(m_mnMesh,n_mnMesh,exy_mn, levels=lvs)
plt.plot(mBaseline,nBaseline,'ko')
plt.clabel(mn_exy_contour, inline=1, fontsize=10)
#plt.title('m-n plane\nShear Strains, Constant a={const}'.format(const=aBaseline))
plt.axis([0,50,0,50])
plt.xlabel('m [# plies]')
plt.ylabel('n [# plies]')
pylab.savefig('mnShear.pdf', bbox_inches='tight')

plt.show()
