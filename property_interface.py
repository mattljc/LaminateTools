"""Superclass for all classes generating material and laminate properties.

This module contains the basic constructors and minimum method definitions
required to define fabric properties and calculate the laminate properties
resulting from this.
"""

from laminate_fundamentals as lf

class Material(object):
    """Superclass for all ply materials.

    Any material (anisotropic, orthotropic, isotropic or otherwise) inherit
    this class.

    Because Hooke's generalized law is used pretty much
    everywhere else in this code they must generate compliance and stiffnes,
    probably as some type of matrix.
    """

    def __init__(self,input_dict):
        # These are the minimum properties that all materials will have.
        self.Name = input_dict['name']
        self.Thickness = input_dict['thk']
        self.Density = input_dict['dens']

    def make_compliance(self):
        """ make_compliance should return a compliance matrix.

        This should be implemented in such a way that is efficient, only
        requiring full calculation the first time the method is called.
        """
        raise NotImplementedError

    def make_stiffness(self):
        """ make_stiffness should return a stiffness matrix.

        Typically the stiffness is the inverse of the compliance. This should
        be implemented in such a way that is efficient, only requiring full
        calculation the first time the method is called.
        """
        raise NotImplementedError

class Properties(object):
    """Superclass for all laminate properties calculations.

    Any property calculator should inherit this class.

    The main goal here is to check that its getting a Laminate type object and
    not just a list of Ply objects as was done in the past. The constructor
    for this also handles total thickness and density. This also defines some
    functions that return useful properties of laminates.
    """

    def __init__(self, laminate_input):
        assert isinstance(lam,lf.Laminate), 'Input not a laminate'
        self.Laminate = laminate_input
        self.TotalThickness = 0
        self.TotalDensity = 0

        for ply in self.Laminate.PlyStack:
            self.TotalThickness += ply.Thickness
            self.TotalDensity += ply.Material.Density

        self.TotalArealDensity = self.TotalDensity * self.TotalThickness

        _ = self.make_global_compliance()
        _ = self.make_global_stiffness()

    def make_global_compliance(self):
        """make_global_compliance should return a compliance matrix that
        converts stress to strain.

        This should be implemented in such a way that is efficient, only
        requiring full calculation the first time the method is called.
        """
        raise NotImplementedError

    def make_global_stiffness(self):
        """make_global_stiffness should return a stiffness matrix that converts
        stress to strain.

        Typically the stiffness is the inverse of the compliance. This should
        be implemented in such a way that is efficient, only requiring full
        calculation the first time the method is called.
        """
        raise NotImplementedError

    def make_effective_properties(self):
        """make_effective_properties should return a dictionary of the
        effective properties for the entire laminate.

        Typically this will be calculated from the elements of the stiffness
        matrix. It may not always be appropriate to do this, for instance in
        asymetric CLPT plates. In such situations warn the user of the problem,
        but do not throw exceptions. Effective properties are very useful for
        sanity checks, and will usually have less than an order of magnitude
        error in them.
        """
        raise NotImplementedError

    def make_strains_from_stress(self, stress):
        """Calculate the strains in laminate resulting from given stresses.

        The stresses should be in the form of a vector. The resulting stresses
        will be in the form of a vector of the same size. This function need
        not be overiden. It should work in the majority of situations, as it
        uses the make_global_stiffness method to calculate the strain.
        """
        self.LaminateStress = stress
        self.LaminateStrain = self.make_global_stiffness() * stress
        return self.LaminateStrain

    def make_stress_from_strains(self, strain):
        """Calculate the stess in laminate resulting from given strains.

        The strains should be in the form of a vector. The resulting stresses
        will be in the form of a vector of the same size. This function need
        not be overiden. It should work in the majority of situations, as it
        uses the make_global_compliance method to calculate the strain.
        """
        self.LaminateStrain = strain
        self.LaminateStress = self.make_global_stiffness() * strain
        return self.LaminateStress

    def make_ply_stress_strain(self):
        """Calculate the stress and strain in each ply along the principal
        fiber directions, returning arrays of ply stresses and strains.

        Sadly this will be different for each type of analysis. The results
        should be stored within each ply.
        """
        raise NotImplementedError

    def make_failure_index(self, type='hoffman'):
        """Calculates a list of failure indices for each ply according to a
        specified index type.

        This should be implemented in such a way that it calculates the failure
        index for a given load condition in each ply and returns the minimum
        failure index for the whole laminate. This function should also be
        capable of calculating the failure index for any arbitrary failure
        type, although Hoffman is used by default.
        """
        raise NotImplementedError

    def __str__(self):
        # Output basic results of the laminate. Ordered keys
        np.set_printoptions(precision=3, linewidth=256, suppress=False)
        output = "== "+self.__class__.__name__+" Properties ==\n"
        keyList = sorted(self.__dict__.keys())
        # Remove keys that aren't useful to report in result string
        while 'specificNT' in keyList: keyList.remove('specificNT')
        while 'laminate' in keyList: keyList.remove('laminate')
        # Build the result string
        for key in keyList:
            output += (key+" = "+str(self.__dict__[key])+"\n")
        return output

    def __repr__(self):
        # Return all the building block representations
        output = '=='+self.__class__.__name__ +'==\n'
        output += repr(self.laminate)
        return output
