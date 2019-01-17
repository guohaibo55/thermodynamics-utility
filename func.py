"""Thermodynamic functions of enthalpy, entropy, specific heat, free energy.

Note the difference between atm and bar: 1 atm = 101.325 kPa; 
1 bar = 100 kPa (by definition). Most published data were referring to
1 atm, although 1 bar is recommended in the JANAF Thermochemical Handbook.
The difference is usually small at logarithm scale: 
log (1 atm / 1 bar)=log(1.01325) ~= 0.01325. 
I am using 1 atm as the standard pressure.   

Created on May 20, 2014
@author: Haibo Guo

Update 2019-1-3: 
* The unit of pressure has changed from bar to MPa.
  This change is to be consistent with the IAPWS's formula.
* The unit of energies (including free energy, enthalpy) has changed from kJ
  to J. Unit conversions have to be done by users.
"""
import numpy as np
from constants import RT, R, P0

def tf_form1(params, T_ref=RT):
    """ Cp(T)=a + b*T + c/T^2 +d*T^2, returns (Cp(T), Delta_H(T), Delta_S(T)).
    
    This is a commonly used function form for calculating isobaric heat 
    capacity. The terms are from $T^{-2}$ to $T^2$ except the $T^{-1}$ term, 
    which appears to be absent in many textbooks. I do not understand why
    the $T^{-1}$ is special and what makes it different from the $T^{-2}$
    term. Other function forms are possible, so I name this one 'form 1'.
    The name 'tf_form1' means thermodynamic function of form 1.
    
    Arguments
        params: an array [a,b,c,d], where a is in the unit of J/mol/K, 
            b in J/mol/K^2, c in J/mol K, and d in J/mol/K^3.
        T_ref: a reference temperature, default is room temperature 
            (298.15 K). Unit: K. 
            Note on T_ref: Why I need a reference temperature?
            Entropy by its definition is zero at the ground state (T=0 K) 
            for pure substances, and does not need a reference state. 
            But, experimental measurements are not feasible for the ground 
            state (T=0 K is not reachable practically), which is outside the 
            valid domain for the Cp function. Therefore, integration from the 
            ground state is only hypothetical. The most convenient practice 
            is to choose a reference temperature, and integrate from this 
            particular temperature. The entropy defined to be zero at this 
            reference temperature. Enthalpy is also defined to be zero at the 
            reference temperature.
    Returns 
        A list of three functions: (Cp, Delta_H, Delta_S), in which
            Cp(T)=a + b*T + c/T^2 +d*T^2. Unit: J/mol/K;
            Delta_H(T, T2=T_ref)= integration of Cp dT from T2 to T. Unit: J/mol.
            Delta_S(T, T2=T_ref)= integration of Cp/T dT from T2 to T. Unit: 
            J/mol/K.
    """
    a, b, c, d = params
    def Cp(T):
        """specific heat at temperature T. Unit: J/mol/K."""
        return a + b*T + c/(T*T) + d * T*T

    def Delta_H(T, T2=T_ref):
        """change of enthalpy from T2 to T. This is integration of Cp dT."""
        return (T*(a+0.5*b*T+1.0/3*d*T*T)-c/T) - (
                T2*(a+0.5*b*T2+1.0/3*d*T2*T2)-c/T2)

    def Delta_S(T, T2=T_ref):
        """change of entropy from T2 to T. This is integration of Cp/T dT."""
        return a*np.log(T/T2) + T*(b+0.5*d*T) - 0.5*c/(T*T) - (
                T2*(b+0.5*d*T2)-0.5*c/(T2*T2))
    # return the three functions in a tuple.
    return (Cp, Delta_H, Delta_S)


def gp_ideal_gas(T, P):
    """ chemical potential's change caused by pressure for ideal gases.

    Arguments
        T: temperature, unit: K.
        P: pressure, unit: MPa.
    Returns
        chemical potential relative to the thermodynamic standard state, 
        according to the equation of state of ideal gases. Unit: kJ/mol.
        The chemical potential is R*T*ln(p/p0).
    """
    return R*T * np.log(P/P0)


class Enthalpy:
    """ Enthalpy of a substance. 
    
    This is a shallow encapsulation of func_H and H_standard.
    The enthalpy at the standard state must be defined.
    The relative enthalpy with reference to the standard state is defined
    in a function.
    It should have functions to
    * calculate relative enthalpy between two temperatures.
    * return the standard-state enthalpy.
    """
    def __init__(self, func_h, std_h):
        """ 
        func_h: a function of temperature, for calculating relative 
             enthalpy; its unit is J/mol.
        std_h: a number; its value is standard-state enthalpy in J/mol.
        """
        self._func_h = func_h
        self._std_h  = std_h

    def h_std(self):
        return self._std_h

    def h_rel(self, T, T_ref=RT):
        """Difference in enthalpy between T_ref to T.
        
        Returns enthalpy different in unit J/mol."""   
        return self._func_h(T, T_ref)
    
    def h(self, T):
        """Standard enthalpy at temperature T.
        
        Returns enthalpy in unit J/mol."""
        return self._func_h(T, RT) + self._std_h

class ChemicalPotential:
    """Chemical potentials of substances.
    
    The central class. The members of a ChemcalPotential object are 
    function of specific heat, function of enthalpy, function of 
    entropy, function of pressure-induced change in chemical potential,
    standard enthalpy of formation and standard entropy.
    
    The use of functions is to make this class applicable to many substances,
    including gases (real or ideal), solids, and liquids, provided the functions
    of the thermodynamic properties can be defined. For example, for solids,
    pressure-induced change in chemical potential can be very small, or 
    negligible, especially when there are also gases in the system.
    
    I need to calculate nonstandard pressure-induced change in chemical
    potential. The assumption of ideal gas should work well for H2, fairly 
    well for O2, but poorly for H2O at high pressures and low temperatures. 
    Laws that work for real gases should be used, instead of those for ideal 
    gases.
    """
    def __init__(self, func_Cp, func_H, func_S, func_gp, params) :
        """ Chemical potential.
        
        Arguments:
        func_Cp: function f(T), heat capacity (in J/mol/K) at T (in K).
        func_H:  function f(T, T_ref), change of enthalpy (J/mol) from T2 
                 (default room temperature, in K) to T (in K).
        func_S:  function f(T, T_ref), change of entropy (J/mol/K) from T2
                 (default room temperature, in K) to T (in K).
        func_gp: function f(T, P), change of chemical potential (in J/mol/K)
                 induced by nonstandard pressure P (in MPa) at a constant 
                 temperature T (in K).
        params: [std_H, std_S], where:
            - std_H is the standard enthalpy (in J/mol), 
            - std_S is the standard entropy (J/mol/K),
        """
        self.func_Cp = func_Cp
        self.func_H = func_H
        self.func_S = func_S
        self.func_gp = func_gp
        self.std_H = params[0]
        self.std_S = params[1]

    def g(self, T, P) :
        """ Free energy of a condition (T,P) relative to the standard state.
        
        Arguments the C API is to thro
            T: temperature (unit: K).
            P: pressure (unit: MPa).
        
        Returns
            Free energy for the condition of T and P, relative to the 
            thermodynamic standard state. Unit: kJ/mol.
        """
        return (self.func_H(T) - T*self.func_S(T) - (T-RT)*self.std_S
                + self.func_gp(T,P)
               )
    

    def gfunc_isobaric(self, P) :
        """ Returns g as a function of temperature at a constant pressure.

        The reference state is the thermodynamic standard state.
        This function returns a function object, which is useful to show 
        temperature dependence of free energy at a constant pressure.
        
        Arguments
            P: pressure, unit MPa.
            
        Returns
            A function object, whose argument is temperature (unit: K), and
            whose return value is free energy (unit: J/mol).
        """
        # define the function to be returned.
        def TD(T) :
            """ returns delta g in J/mol."""
            return (self.func_H(T) - T*self.func_S(T) - (T-RT)*self.std_S
                     + self.func_gp(T,P)
                   )
        return TD

    
    def gfunc_isothermal(self, T) :
        """ Returns g as a function of pressure at a constant temperature.

        The reference state is the thermodynamic standard state.
        This function returns a function object, which is useful to show 
        pressure dependence of free energy at a constant temperature.
        
        Arguments
            T: temperature, unit K.
            
        Returns
            A function object, whose argument is pressure (unit: MPa), and
            whose return value is free energy (unit: J/mol).
        """
        term_T = self.func_H(T) - T*self.func_S(T) - (T-RT)*self.std_S
        # defines the function
        def PD(P) :
            """ returns delta mu in J/mol."""
            return term_T + self.func_gp(T, P)
        return PD

    def iso_mu_IG(self, m):
        """ returns a function P(T) corresponds to a chemical potential.

        m: chemical potential in J/mol.
        Returns:
        PP(T): partial pressure as a function of temperature.
            PP in bar,
            T in K.
        """
        def PP(T):
            term_T = self.func_H(T) - T*self.func_S(T) - (T-RT) * self.std_S
            term_P = np.exp((m - term_T) / (R*T))
            return term_P
        return PP

# Testing
if __name__ == "__main__" :
    # parameters, be cautious about units, valid range
    _tf_H2  = tf_form1([31.31895, -0.005145276, -119747.5,  4.102189e-6])
    _tf_O2  = tf_form1([23.93617,  0.01650568,   86891.96, -5.626105e-6])
    _tf_H2O = tf_form1([29.34480,  0.01086696,   73907.32,  9.060696e-7])
    
#    pot_H2 = ChemicalPotential(_tf_H2[0], _tf_H2[1], _tf_H2[2], gp_ideal_gas, 
#                     [0.0, 130.68, -695.602])
#    pot_O2 = ChemicalPotential(_tf_O2[0], _tf_O2[1], _tf_O2[2], gp_ideal_gas, 
#                     [0.0, 205.147, -1003.227])
#    pot_H2O = ChemicalPotential(_tf_H2O[0], _tf_H2O[1], _tf_H2O[2], gp_ideal_gas, 
#                     [0.0, 188.834, -1417.455])

