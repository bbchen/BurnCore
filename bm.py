#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Birch Murnaghan EOS module
Author: Bin Chen
"""

def bm_eos(k0, k0p, v0_v):
    """Birch-Murnaghan equation of stae
    P(V) = 1.5 *k0 * (DensityRatio^(7.0/3.0)- (DensityRatio(^5.0/3.0)) * (1+0.75 * (k0P-4) * (Density^(2/3)-1))

    Parameters:

    :param k0: bulk modulus at zero pressure
    :param k0p: pressure derivative of bulk modulus
    :param v0_v: ratio of zero-pressure volume over volume
    :returns: pressure at $V_0/V$

    """
    Pressure = 1.5*k0*(v0_v**(7.0/3.0)-v0_v**(5.0/3.0)) * \
          (1 + 0.75*(k0p-4)*(v0_v**(2.0/3.0)-1))
    return Pressure

def v_bm(k0, k0p, v0, pressure):
    """Calculate V from the Birch-Murnaghan EOS."""
    v0_v = (1.0  + pressure*k0p/k0)**(1.0/k0p)    #v0_v = V0/V

    # iteration
    tol = 1e-9
    nmax = 1e4
    n = 0
    pressure_try = bm_eos(k0, k0p, v0_v)
    while ((n<nmax) and abs(pressure_try - pressure) > tol ):
        v0_v += (pressure - pressure_try)/(k0+k0p*pressure)*v0_v
        n += 1
        if (n == nmax):
          v0_v = 0
        pressure_try = bm_eos(k0, k0p, v0_v)
    # Calc Volume
    volume = v0 / v0_v
    return volume

def k_bm(k0, k0p, v0_v):
    """Calculate the B(P)

    :param v0_v: V0/V

    """
    Ta = 0.5 *k0*( 7 * v0_v**(7.0/3.0)- 5.0 * v0_v**(5.0/3.0))
    Tb = 3 /8.0 *k0*(k0p - 4)*(9.0 * v0_v**3.0 - 14 * v0_v**(7.0/3.0) + 5.0*v0_v**(5.0/3.0))
    B = Ta + Tb
    return B

def kp_bm(k0, k0p, v0_v):
    """
    Calculate the B'(P) (ref : Singh2005a)
    v0_v : V0/V
    """
    kp = k0/(8.0*k_bm(k0,k0p,v0_v)) * \
          ((k0p-4.0)*(81 * v0_v**3.0  - 98 * v0_v**(7.0/3.0) \
          + 25.0*v0_v**(5.0/3.0)) + 4.0/3.0 *(49.0*v0_v**(7.0/3.0 ) \
          - 25.0 *v0_v**(5.0/3.0 )))
    return kp