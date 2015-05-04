import numpy as np
import scipy.optimize as opti
import matplotlib
import matplotlib.cm as cm
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

import materials
import lamination
import constants
import plate_failure
import failure_envelope


# Material Definition
E1 = 21.3e6
E2 = 1.5e6
Nu12 = 0.27
G12 = 1.0e6
CPT = 0.01
F1T = 330e3
F1C = -250e3
F2T = 8.3e3
F2C = -33e3
F12S = 11e3

propDict = {'name':'uni','Dens':1.0,'CPT':CPT,'E11':E1,'E22':E2,
            'Nu12':Nu12,'G12':G12,'f1t':F1T,'f1c':F1C,'f12s':F12S,
            'f2t':F2T,'f2c':F2C}
uni = materials.Plate(propDict)
core = materials.WeakCore()

# Validation run
plyValid = lamination.Ply({'matl':uni,'thk':0.01,'orient':0.0})
lamValid = lamination.Laminate(plyBook=[plyValid], symmetry=True, n_count=2)
propsValid = constants.ThinPlate(lamValid)
#print(propsValid.A.I)
#print(plyValid.Material.Compliance.I)

loadsValid = np.matrix([1000,0,0,0,0,0]).T
propsValid.calculatePlyStressStrain(loadsValid)
stressValid = plate_failure.MaxStress(propsValid.laminate)
hoffmanValid = plate_failure.Hoffman(propsValid.laminate)
#print(stressValid)
#print(hoffmanValid)


# Define laminates and analysis objects
ply0thin = lamination.Ply({'matl':uni,'thk':0.01,'orient':0.0})
ply30thin = lamination.Ply({'matl':uni,'thk':0.01,'orient':30.0})
plyN30thin = lamination.Ply({'matl':uni,'thk':0.01,'orient':-30.0})
ply90thin = lamination.Ply({'matl':uni,'thk':0.01,'orient':90.0})
plybookThin = [ply0thin, ply30thin, plyN30thin, ply90thin]

ply0sandwich = lamination.Ply({'matl':uni,'thk':0.01,'orient':0.0})
ply30sandwich = lamination.Ply({'matl':uni,'thk':0.01,'orient':30.0})
plyN30sandwich = lamination.Ply({'matl':uni,'thk':0.01,'orient':-30.0})
ply90sandwich = lamination.Ply({'matl':uni,'thk':0.01,'orient':90.0})
plyCore = lamination.Ply({'matl':core,'thk':0.25,'orient':0.0})
plybookSandwich = [ply0sandwich, ply30sandwich, plyN30sandwich, ply90sandwich, plyCore]

thin = lamination.Laminate(plyBook=plybookThin, symmetry=True)
sandwich = lamination.Laminate(plyBook=plybookSandwich, symmetry=True)

thinProps = constants.ThinPlate(thin)
sandwichProps = constants.ThinPlate(sandwich)

#ResThin = makeEnvelope(thinProps,plate_failure.Hoffman,'any')
ResSand = failure_envelope.makeLinVarEnvelope(sandwichProps,plate_failure.Hoffman,'any')

aPlot = [ResSand['aMin']] + ResSand['aaRange'] + [ResSand['aMax']] + ResSand['aaRange'][::-1] + [ResSand['aMin']]
bPlot = [0] + ResSand['bUpper'] + [0] + ResSand['bLower'][::-1] + [0]

plt.figure(1)
plt.plot(aPlot,bPlot)
plt.axis([-2,2,-4,4])
plt.grid(b=True, which='major', color='k', linestyle='--')
plt.xlabel('Force Resultant Multiplier, a')
plt.ylabel('Moment Resultant Multiplier, b')
plt.show()
