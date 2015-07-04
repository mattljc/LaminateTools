import warnings
import numpy as np
import property_interface

class Plate(property_interface.Material):
    """A plate material for use in classical laminated plate theory (CLPT).

    In addition to the terms required for the Material superclass, CLPT
    Plate materials require the 4 in-plane properties: E11, E22, Nu12 and G12.
    See wiki for more details on theoretical aspects.
    """

    def __init__(self, property_dict = None):
        property_interface.Material.__init__(self, property_dict)

        # MINIMUM FUNCTIONAL PROPERTY SET
        try:
            self.E11 = float(property_dict['E11'])
            self.E22 = float(property_dict['E22'])
            self.Nu12 = float(property_dict['Nu12'])
            self.G12 = float(property_dict['G12'])
        except KeyError:
            raise KeyError('Check input, minimum information not provided')

        # These properties are not required for basic functionality.
        # CTE
        try:
            self.A1 = propsDict['a1']
            self.A2 = propsDict['a2']
        except KeyError:
            warnings.warn('No CTE included, setting all CTE to zero')
            self.A1 = 0
            self.A2 = 0

        # STRESS LIMITS
        try:
            self.F1t = float(self.InputDict['f1t'])
            self.F1c = float(self.InputDict['f1c'])
            self.F2t = float(self.InputDict['f2t'])
            self.F2c = float(self.InputDict['f2c'])
            self.F12s = float(self.InputDict['f12s'])
        except KeyError:
            warnings.warn('No stress limits included, setting all to infinity')
            self.F1t = np.inf
            self.F1c = np.inf
            self.F2t = np.inf
            self.F2c = np.inf
            self.F12s = np.inf

        # STRAIN LIMITS
        try:
            self.Ep1t = self.InputDict['e1t']
            self.Ep1c = self.InputDict['e1c']
            self.Ep2t = self.InputDict['e2t']
            self.Ep2c = self.InputDict['e2c']
            self.Ep12s = self.InputDict['e12s']
        except KeyError:
            warnings.warn('No strain limits included, setting all to infinity')
            self.Ep1t = np.inf
            self.Ep1c = np.inf
            self.Ep2t = np.inf
            self.Ep2c = np.inf
            self.Ep12s = np.inf

        # TODO: Add necessary handlers for vibration analysis

    def make_compliance(self):
        try:
            return self.Compliance
        except AttributeError:
            s11 = 1/self.E11
            s12 = -self.Nu12/self.E11
            s22 = 1/self.E22
            s66 = 1/self.G12
            self.Compliance = np.matrix([ \
            [s11, s12, 0  ], \
            [s12, s22, 0  ], \
            [0  , 0  , s66]])
            return self.Compliance

    def make_stiffness(self):
            self.Stiffness = self.make_compliance.I

    def make_invariants(self):
        try:
            return [self.U1, self.U2, self.U3, self.U4, self.U5]
        except AttributeError:
            self.U1 = (Q[0,0] + Q[1,1])*3/8 + Q[0,1]/4 + Q[2,2]/2
            self.U2 = (Q[0,0] - Q[1,1])/2
            self.U3 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 - Q[2,2]/2
            self.U4 = (Q[0,0] + Q[1,1])/8 + Q[0,1]*3/4 - Q[2,2]/2
            self.U5 = (Q[0,0] + Q[1,1])/8 - Q[0,1]/4 + Q[2,2]/2
            return [self.U1, self.U2, self.U3, self.U4, self.U5]

class ThinPlates(property_interface.Properties):
    """Calculates the elastic properties, CTE and dynamic response of a thin
    laminate where out-of-plane properties can be ignored using Classic
    Laminated Plate Theory (CLPT).

    The assumption of infinitessimal transverse shears requires that the ratio
    of thickness to the short transverse length be less than 1:10. If this is
    not true for the laminate under consideration.

    Note that z-coordinates are zero at the mid-surface as expected from CLPT
    and are negative in the direction towards the tool surface.
    """

    def __init__(self, lam=None):
        try:
            property_interface.Properties.__init__(self,lam)
        except AttributeError:
            raise  AttributeError('Check material has the required properties')

    def make_global_stiffness(self):
        """Builds ABD matrix ply by ply, and returns the augmented ABD matrix.

        In CLPT the concept of a stiffness matrix is hard to define by itself.
        The closest representation is the A matrix, but is not useful in cases
        where there is substantial bending. Instead, the full ABD matrix is
        returned to the user.

        This function also generates thermal and dynamic properties.
        """
        try:
            return self.ABD
        except AttributeError:
            self.A = np.matrix( np.zeros((3,3)) )
            self.B = np.matrix( np.zeros((3,3)) )
            self.D = np.matrix( np.zeros((3,3)) )
            self.ABD = np.matrix( np.zeros((6,6)) )
            self.specificNT = np.matrix( np.zeros((3,1)) )

            zLow = -self.TotalThickness / 2

            # Build global stiffness ply by ply, then add to the ABD matrix.
            # Invariant method is used to ensure matrix symetry in the final result.
            for ply in self.laminate.PlyStack:
                zUp = zLow + ply.Thickness
                U = ply.Material.make_invariants()
                c4 = np.cos(4 * np.radians(ply.Orientation))
                c2 = np.cos(2 * np.radians(ply.Orientation))
                c1 = np.cos(1 * np.radians(ply.Orientation))
                s4 = np.sin(4 * np.radians(ply.Orientation))
                s2 = np.sin(2 * np.radians(ply.Orientation))
                s1 = np.sin(1 * np.radians(ply.Orientation))

                Q11 = U[0] + U[1]*c2 + U[2]*c4
                Q22 = U[0] - U[1]*c2 + U[2]*c4
                Q12 = U[3] - U[2]*c4
                Q66 = U[4] - U[2]*c4
                Q16 = U[1]*s2/2 + U[2]*s4
                Q26 = U[1]*s2/2 - U[2]*s4
                ax = ply.Material.a1*c1**2 + ply.Material.a2*s1**2
                ay = ply.Material.a1*s1**2 + ply.Material.a2*c1**2
                axy = (ply.Material.a2-ply.Material.a1)*c1*s1

                ply.GlobalStiffness = np.matrix([
                [Q11, Q12, Q16], \
                [Q12, Q22, Q26], \
                [Q16, Q26, Q66]])
                ply.GlobalCTE = np.matrix([[ax],[ay],[axy]])

                self.A += ply.GlobalStiffness * (zUp - zLow)
                self.B += ply.GlobalStiffness * (zUp**2 - zLow**2) / 2
                self.D += ply.GlobalStiffness * (zUp**3 - zLow**3) / 3
                self.specificNT += ply.GlobalStiffness * ply.GlobalCTE * (zUp - zLow)

                # Increment Z
                zLow = zUp

            self.ABD[0:3,0:3] = self.A
            self.ABD[0:3,3:] = self.B
            self.ABD[3:,0:3] = self.B
            self.ABD[3:,3:] = self.D

    def make_global_compliance(self):
        """Returns the inverted ABD matrix.

        This is also a tricky one to define for CLPT laminates. For similar
        reasons to make_global_stiffness this returns the inveted ABD matrix
        instead of just the global compliance matrix.
        """
        return self.make_global_stiffness.I

    def make_effective_properties(self):
        """Returns a dictionary containing overall laminate properties. This
        dictionary will contain thermal and dynamic properties also.

        Note that for CLPT, these properties are only valid for symetric or
        pseudo-symetric laminates. Anything else will generate a warning.
        """
        try:
            return self.EffectiveProperties
        except AttributeError:
            effective_compliance = self.make_global_compliance()
            Exx = 1 / (effective_compliance[0,0] * self.TotalThickness)
            Eyy = 1 / (effective_compliance[1,1] * self.TotalThickness)
            Gxy = 1 / (effective_compliance[2,2] * self.TotalThickness)
            Nuxy = - effective_compliance[0,1] / effective_compliance[0,0]
            Etaxs = effective_compliance[0,2] / effective_compliance[0,0]
            Etays = effective_compliance[1,2] / effective_compliance[1,1]

            effective_CTE = effective_compliance * self.specificNT
            ax = effective_CTE[0,0]
            ay = effective_CTE[1,0]
            axy = effective_CTE[2,0]

            self.EffectiveProperties = {'Exx':Exx,
                                        'Eyy':Eyy,
                                        'Gxy':Gxy,
                                        'Nuxy':Nuxy,
                                        'Etaxs':Etaxs,
                                        'Etays':Etays,
                                        'ax':ax,
                                        'ay':ay,
                                        'axy':axy}
            return self.EffectiveProperties

    def make_strains_from_stress(self, resultants):
        """Calculate global mid-plane strains and curvatures given force and
        moment resultants (through thickness integrated stress) augmented as a
        single vector.
        """
        self.Resultants = resultants
        self.StrainsCurves = self.ABD.I * self.Resultants
        return self.StrainsCurves

    def make_stress_from_strains(self, strains_curves):
        """Calculate global thickness integrated stresses (resultants) from
        global midsurface strains and curvatures.
        """
        self.StrainsCurves = strains_curves
        self.Resultants = self.ABD * self.StrainsCurves
        return self.Resultants

    def make_ply_stress_strain(self, resultants=None, strain=None):
        """Using the global strains, determines the stress and strains in each
        ply of the laminate.

        This function returns arrays of stress and strain in fibre coordinates
        on a ply by ply basis.
        """
        try:
            midStrains = self.StrainsCurves[:3,0]
            midCurves = self.StrainsCurves[3:,0]
        except AttributeError:
            if not strain==None:
                self.make_stress_from_strains(strain)
            elif not stress==None:
                self.make_strains_from_stress(resultants)
            midStrains = self.StrainsCurves[:3,0]
            midCurves = self.StrainsCurves[3:,0]

        stress_array = np.zeros((len(self.Laminate.PlyStack),3))
        strain_array = np.zeros((len(self.Laminate.PlyStack),3))
        counter = 0

        zLow = -self.TotalThickness / 2
        for ply in self.Laminate.PlyStack:
            m = np.cos(np.radians(ply.Orientation))
            n = np.sin(np.radians(ply.Orientation))
            T = np.matrix([\
            [m**2, n**2, 2*m*n],\
            [n**2, m**2, -2*m*n],\
            [-m*n, m*n, m**2-n**2]])

            zUp = zLow + ply.Thickness
            # Uses ply mid-surface to determine state of each ply.
            globalStrain = midStrains + (zUp+zLow)/2*midCurves
            globalStrain[2,0] = globalStrain[2,0]/2 #Convert gamma to epsilon

            ply.Strain = T * globalStrain
            ply.Stress = ply.Material.Compliance.I * ply.Strain
            zLow=zUp

            stress_array[counter,:] = ply.Stress
            strain_array[counter,:] = ply.Strain

        return (stress_array, strain_array)

    def make_failure_index(self, type='hoffman'):
        index = list()

        for ply in self.Laminate.PlyStack:
            s1 = ply.Stress[0,0]
            s2 = ply.Stress[1,0]
            s12 = ply.Stress[2,0]
            eps1 = ply.Strain[0,0]
            eps2 = ply.Strain[1,0]
            eps12 = ply.Strain[2,0]

            if type=='hoffman':
                f1t = ply.Material.f1t
                f1c = ply.Material.f1c
                f2t = ply.Material.f2t
                f2c = ply.Material.f2c
                f12s = ply.Material.f12s
                ply.HoffmanFail = -s1**2/(f1t*f1c) + s1*s2/(f1t*f1c) \
                                -s2**2/(f2t*f2c) + s1*(1/f1t+1/f1c) \
                                +s2*(1/f2t+1/f2c) + (s12/f12s)**2
                index.append(ply.HoffmanFail)

            elif type=='tsaihill':
                f1 = ply.Material.f1t*(s1>=0) + ply.Material.f1c*(s1<0)
                f2 = ply.Material.f2t*(s2>=0) + ply.Material.f2c*(s2<0)
                f12 = ply.Material.f12s
                ply.TsaiHillFail = (s1/f1)**2 + (s2/f2)**2 + (s12/f12)**2 \
                                   -s1*s2/(f1**2)
                index.append(ply.TsaiHillFail)

            elif type=='maxstress':
                f1 = ply.Material.f1t*(s1>=0) + ply.Material.f1c*(s1<0)
                f2 = ply.Material.f2t*(s2>=0) + ply.Material.f2c*(s2<0)
                f12 = ply.Material.f12s
                ply.MaxStressFail = [s1/f1, s2/f2, s12/f12]
                index.append(ply.MaxStressFail)

            elif type=='maxstrain':
                e1 = ply.Material.e1t*(e1>=0) + ply.Material.e1c*(e1<0)
                e2 = ply.Material.e2t*(e2>=0) + ply.Material.e2c*(e2<0)
                e12 = ply.Material.e12s
                ply.MaxStrainFail = [eps1/e1, eps2/e2, eps12/e12]
                index.append(ply.MaxStrainFail)

        return np.array(index)

if __name__=="__main__":
    matl_dict = {'name':'AS4-8552-UNI',
                'thk':0.0074,
                'dens':0.057,
                'E11':19.09e6,
                'E22':1.34e6,
                'Nu12':0.335,
                'G12':0.70e6,
                'f1t':279.61e3,
                'f1c':215.29e3,
                'f2t':9.27e3,
                'f2c':38.85e3,
                'f12s':13.28e3,}
    matl = thin_plates.Plate(matl_dict)

    plate00 = laminate_fundamentals.Ply({'matl':matl,'thk':0.0074,'orient':0})
    plate30 = laminate_fundamentals.Ply({'matl':matl,'thk':0.0074,'orient':30})
    plate60 = laminate_fundamentals.Ply({'matl':matl,'thk':0.0074,'orient':60})
    plate90 = laminate_fundamentals.Ply({'matl':matl,'thk':0.0074,'orient':90})
    stack = [plate00,plate30,plate60,plate90]

    lam = laminate_fundamentals.Laminate(stack, 2, True)
    res = thin_plates.ThinPlates(lam)
