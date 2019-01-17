'''define physical constants, unit converting factors.
Created on Jun 19, 2011

@author: haibo
'''

Boltzmann = 1.38064852E-23   # Boltzmann constant, J/K.
Planck = 6.626070040E-34     # Planck constant, J S.
Planck_reduced = 1.054571800E-34   # reduced Planck constant, J S.

# Atomic unit
AMU = 1.660539040E-27     # atomic mass unit, kg.

# Thermodynamics
RT = 298.15    # room temperature, K
R = 8.314      # gas constant, J/mol/K

N0 = 6.0221415E23    # avogadro number, dimensionless.
N0_base = 6.0221415    # avogadro number without exponential

P0 = 0.101325    # Mpa

eV_kJ = 94.4861489     # eV per atom to kJ per mol.
eVperA2 = 16.0219     # convert the surface energy in eV/A2 to J/m2

SVP_H2O_RT = 3.2/101.325   # saturated vapor pressure of water at RT
