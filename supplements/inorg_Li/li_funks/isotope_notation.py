'''
Functions for converting between isotope notation.

(c) 2019 : Oscar Branson : oscarbranson@gmail.com
'''

def R_2_A(R):
    return R / (1 + R)

def A_2_R(A):
    return A / (1 - A)

def R_2_alpha(Rsample, Rstd):
    return Rsample / Rstd

def alpha_2_R(alpha, Rstd):
    return alpha * Rstd

def alpha_2_A(alpha, Rstd):
    return R_2_A(alpha_2_R(alpha, Rstd))

def A_2_alpha(A, Rstd):
    return R_2_alpha(A_2_R(A), Rstd)