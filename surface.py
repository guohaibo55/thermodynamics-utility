''' Thermodynamics of surfaces.
A surface may have multiple terminations, each having its own stoichiometry, 
surface energy, and surface stress.

Created on June 7, 2018

@author: Haibo Guo (guohaibo@shu.edu.cn)
'''
import numpy as np
from operator import itemgetter
from copy import copy

import constants

class Surface:
    """ A surface that may have multiple terminations. """
    def __init__(self, area, label="", conversion=0.166053886313E-3):
        """ Define properties: area, label.

        Arguments
        area: (scalar; float)
            Surface unit cell's area. scalar, float. Its unit is not
            hard-coded at this time; but it must be consistent with surface
            energies. Candidate for its unit is: A2 (angstrom squared), m2, or cm2.
        label: (scalar; string)
            User-defined label for the surface, e.g. "(100)",
            "(1 0 -1 0)", "Fe3O4(111)", etc.
        conversion: (scalar; float)
            This is to convert the unit of energy/area to J/m2.
            The default value, 0.166053886313E-3, is for the case when
            delta_mu is in J/mol, and area is in A2.                
        """
        # I should check that area>0.0 here.
        self.area = area
        self.label = copy(label)
        self.terminations = []
        self.conversion = conversion

    def push_termination(self, energy, excess, label="", ext=[]):
        """ New termination.

        This function is to append a termination to the terminations list.
        There is no return value.
        Arguments
            energy: surface energy for this termination at a reference state.
                Its unit is not hard-coded, but should be consistent with the
                Surface object's area. Frequently used units are: eV/A2, J/m2.
                (float, scalar)
            excess: numbers of excess species per unit-cell. Dimensionless.
                The numbers may not be integers. For a stoichiometric
                termination, all elements of the excess list is 0.
                (float, list)
            label: name of the termination.
                (string, scalar)
        """
        self.terminations.append(
            Termination(energy, self.area, excess, label, self.conversion, ext)
            )
                

    def min_termination(self, delta_mu):
        """ Calculates the minimum-energy terminations for a given conditions.

        Arguments
        delta_mu: (list; float)
            a list of chemical potentials of excess spcies.
            The unit must be compatible with the unit conversion factor.

        Returns
        [gamma, index]: the minimum surface energy among the terminations,
            and the index corresponding to minimum-energy termination for
            the chemical potential array delta_mu.

        Raises
        ValueError if len(delta_mu) does not match that of self.excess.
        """
        if len(self.terminations)==0 :
            print("There is no terminations. Returning [None, None].")
            return [None, None]
        else :
            if len(self.terminations[0].excess) != len(delta_mu):
                print "Error: incompatible dimension of excess and delta_mu"
                raise ValueError
            gammas = [t.gamma(delta_mu) for t in self.terminations]
            idx_min = min(enumerate(gammas), key=itemgetter(1))[0]
            return [gammas[idx_min], idx_min]


class Termination :
    """ A surface termination. 
    
    Data required to construct a surface termination is stoichiometry, energy, 
    area, label, and extra data.
    """
    def __init__(self, energy, area, excess,
                 label="", conversion=0.166053886313E-3, ext=[]) :
        """ Define the termination's properties. 

        Arguments
        energy: (scalar; float).
            Surface excess energy at thermodynamic standard state (unit: energy
            per area, such as J/m2 or eV/A2; unsettled).
            This energy should be normalized per area.
        area: (scalar; float) 
            Area of a surface unit cell (unit of length^2; unsettled). 
            Tis area is to divide chemical potential, thus converting 
            energy (unit: J, cal, eV, etc.) to energy per area (J/m2, eV/A2,
            etc.). 
        excess: (list; float)
            numbers of excess species (mostly molecules) per unit-cell
            (dimensionless). This is a list. Multiple types of molecules can be
            put into this list. If there is only one type of molecules, it
            should be put into a list anyway.
        label: (scalar; string)
            a brief description (name, label) of the termination.
        conversion: (scalar; float)
            This is to convert the unit of energy/area to J/m2 when delta_mu is
            in J/mol, and area is in A2.            
        ext: (list; undefined)
            a list of extended properties of the surface. 
            The interpretation of the properties are to be defined by user. 
            This allows user-defined function to specifically interpret the 
            properties. For example, one may add number of surface hydroxyl 
            groups to this list.

        Notes on unsettled units and unit conversion:
            I do not hard code the units for the argumenst energy and area.
            The calling program should do the job of unit conversion. If
            area is in angstrom^2 and delata_mu is in J/mol, then the number
            delta_mu/area should be multiplied by conversion to convert its
            unit into J/m2.
        """
        self.energy = energy
        self.area = area
        self.excess = np.array(excess)
        self.label = copy(label)
        self.conversion = conversion
        self.ext = copy(ext)
    
    def gamma(self, delta_mu):
        """ Calculates surface energy at given chemical potentials.

        Arguments
        delta_mu: (list; float)
            an array of chemical potential (in J/mol) with reference
            to the thermodynamic standard state. Its size (or length) must
            be equal to that of excess.

        Returns
            Surface energy at new chemical potentials (in J/m2).

        Raises
        ValueError if len(delta_mu) does not match that of self.excess.
        """
        if len(delta_mu) != len(self.excess) :
            print "Incompatible numbers of excess species and delta_mu."
            raise ValueError
        delta_gamma = 0.0
        for i,t in enumerate(self.excess):
            delta_gamma += t*delta_mu[i]
        delta_gamma = delta_gamma/(2.0*self.area) * self.conversion
        return self.energy - delta_gamma
    

def minE_termination(terminations, delta_mu):
    """ Calculates the minimum-energy terminations under the given conditions.
    
    terminations: a group of terminations being compared with one another.
                  It is not very meaningful to compare different surface 
                  orientations, such as comparing hematite (012) and (001) 
                  surfaces. Instead, it is meaningful to compare different 
                  terminations of a same orientation, such as to find the 
                  minimum-energy termination of hematite (012) under a given
                  condition.
    Returns: [gamma, index] corresponding to the condition.
    """
    if len(terminations)<2 : 
        print "Warning: there is only one or zero terminations.", \
             "No comparisons are necessary."
    if len(terminations[0].excess) != len(delta_mu):
        print "Error: incompatible dimension of excess and delta_mu"
        raise ValueError

    gammas = [t.gamma(delta_mu) for t in terminations]
    idx_min = min(enumerate(gammas), key=itemgetter(1))[0]
    return [gammas[idx_min], idx_min]

