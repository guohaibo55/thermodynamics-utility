""" Define a set of phases for which thermodynamic functions are defined.
Among the thermodynamic functions, free energy or chemical potential is
essential; other functions, such as enthalpy, entropy, and heat capacity
are to be added. It's difficult to invent names for these phases that are
consistent, and here is a tentative scheme: formula_state, where state is
one of G (for gas), S (for solid), L (liquid), LG (liquid and gas). 
For example, O2_G, H2O_LG."""

from iapws import IAPWS95

import func


class phase_base:
    """ Base class of phases.

    This class contains abstract methods that are to be implemented by
    subclasses. It only defines interfaces (function names, arguments, 
    units of input and output arguments), which should be abide by
    subclasses. At present (2018-1-7) only g(T,P) is mandatory to 
    be implemented, but more thermodynamic functions will be included.
    """

    def __init__(self):
        pass

    def g(self, T, P):
        """ Calculate free energy for a temperature and pressure.

        Arguments
            T: temperature (unit: K);
            P: pressure (unit: MPa);
        Returns
            Free energy in J/mol.
        """
        pass


class O2_G(phase_base):
    """ The O2 gas phase in temperature range 250~1200 K.
    """
    def __init__(self):
        tf  = func.tf_form1([23.93617, 0.01650568, 86891.96, -5.626105e-6])
        self.g_func = func.ChemicalPotential(
            tf[0], tf[1], tf[2], 
            func.gp_ideal_gas, 
            [0.0, 205.147]
        )

    def g(self, T, P):
        """ Free energy at T and P."""
        return self.g_func.g(T, P)
    
        
    
class H2_G(phase_base):
    """ The H2 gas phase in temperature range 250~1200 K. """
    def __init__(self):
        tf = func.tf_form1([31.31895, -0.005145276, -119747.5,  4.102189e-6])
        self.g_func = func.ChemicalPotential(
            tf[0], tf[1], tf[2], 
            func.gp_ideal_gas, 
            [0.0, 130.68]
        )

    def g(self, T, P):
        """ Free energy at T and P."""
        return self.g_func.g(T, P)
    


class H2O_ideal_gas_G(phase_base):
    """ H2O gas (water vapor) phase in temperature range 200~1200 K. 

    Note: The gas phase of H2O is thermodynamically unstable for temperatures
    below 100 degree celsius under standard pressure. The treatment here is
    to approximate H2O gas as an ideal gas. I have checked the error caused
    by this unrealistic approximation, and found it is big."""
    def __init__(self):
        tf = func.tf_form1([31.31895, -0.005145276, -119747.5,  4.102189e-6])
        self.g_func = func.ChemicalPotential(
            tf[0], tf[1], tf[2], 
            func.gp_ideal_gas, 
            [-285.830*1000, 69.950]
        )

    def g(self, T, P):
        """ Free energy at T and P."""
        return self.g_func.g(T, P)


class H2O_IAPWS95_LG(phase_base):
    """ Liquid and gas phases of H2O defined in IAPWS95 document.
    
    This method requires a module named iapws, that implement IAPWS
    equations. The equation set IAPWS95 is for scientific usage, and
    the equation set IAPWS97 is for industrial usage. The IAPWS95 
    equations rely on IAPWS97, therefore the two equation sets share
    the same domain: for T in [273.15, 1073.15], P<100 MPa; and for
    T in [1073.15, 
    """
    def __init__(self):
        self.g0 = -82.088162723316 #This is IAPWS95(T=298.15, P=0.101325).g*18

    def g(self, T, P):
        # 18.0 is to convert kJ/kg to J/mol
        return IAPWS95(T=T, P=P).g * 18.0 - self.g0
        

class F2_G(phase_base):
    """ The F2 gas phase in temperature range 250~1200 K. """
    def __init__(self):
        tf = func.tf_form1([29.82499449, 0.0102094987, -142876.83,  4.6093085e-6])
        self.g_func = func.ChemicalPotential(
            tf[0], tf[1], tf[2], 
            func.gp_ideal_gas, 
            [0.0, 31.302]
        )

    def g(self, T, P):
        """ Free energy at T and P."""
        return self.g_func.g(T, P)


#=== testing
if __name__ == "__main__" :
    import matplotlib.pyplot as plt
    import numpy as np

    test_P = np.array([100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10])
    test_T = np.array([273.15, 298.15, 373.15, 473.15, 573.15, 673.15])
    
    print "A few chemical potentials of H2, O2, and H2O(gas)."
    O2 = O2_G()
    H2 = H2_G()
    H2O = H2O_IAPWS95_LG()
    
    print "H2: (unit J/mol)\nT=", test_T
    for p in test_P : 
        for x in test_T:
            print "T=", x, "P=", p, "g=", H2.g(x,p)
    
    print "O2:"
    for p in test_P : 
        for x in test_T:
            print "T=", x, "P=", p, "g=", O2.g(x,p)
            
    print "H2O(gas)"
    #for p in test_P :
    #    for x in test_T:
    #        print "T=", x, "P=", p, "g=", H2O.g(x,p)
    temp = np.arange(273.15, 1273.15, 10)
    p1 = 0.0031    # saturated vapor pressure at 298.15 K, 3.1 KPa.
    p2 = 20.0      # typical pressure in hydrothermal autoclave
    a1 = []
    a2 = []
    for t in temp:
        a1.append(H2O.g(t, p1))
        a2.append(H2O.g(t, p2))

    plt.plot(temp, a1, label="SVP")
    plt.plot(temp, a2, label="HT")
    plt.xlabel("Temperature (K)")
    plt.ylabel("G (J/mol)")
    plt.legend()
    plt.show()


