'''
Functions for calculating the d13C of aqueous carbon species and calcite.

(c) 2019 : Oscar Branson : oscarbranson@gmail.com
'''

import uncertainties as un

###########################################
# Aqueous Carbon Species
###########################################

# From Zhang et al. (1995; https://doi.org/10.1016/0016-7037(95)91550-D)

def epsilon_CG(T_C):
    """
    Calculate epsilon_CO3-CO2(g).
    """
    m = un.ufloat(-0.052, 0.021)
    c = un.ufloat(7.22, 0.38)
    return m * T_C + c

def alpha_CG(T_C):
    """
    Calculate alpha_CO3-CO2(g).
    """
    return 1 + epsilon_CG(T_C) / 1000

def epsilon_BG(T_C):
    """
    Calculate epsilon_HCO3-CO2(g).
    """
    m = un.ufloat(-0.1141, 0.0028)
    c = un.ufloat(10.78, 0.04)
    return m * T_C + c

def alpha_BG(T_C):
    """
    Calculate alpha_HCO3-CO2(g).
    """
    return 1 + epsilon_BG(T_C) / 1000

def alpha_BC(T_C):
    """
    Calculate alpha_HCO3-CO3.
    """
    return alpha_BG(T_C) / alpha_CG(T_C)

###########################################
# Mineral Species
###########################################

def alpha_SB():
    """
    Calculate alpha_solid-HCO3.
    """
    return 1 + un.ufloat(1, 0.2) / 1000