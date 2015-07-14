# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:45:57 2013

@author: binchen
"""
# module mgd
""" Mie-Gruneisen Debye Equation of State"""
from math import *
from eos import *

R = 0.00831447 # cm^3 GPa/K/mole

def debye_t(debye_t0, gamma0, gamma_p, q):
    """Calculate Debye temperature."""
    debye_t = float(debye_t0)*exp((float(gamma0) - float(gamma_p)) / float(q))
    return debye_t

def gamma_p(gamma0, q, v0v_ratio):
    """Calculate gamma at high pressure"""
    # v0v_ratio = V0 / V
    gamma_p = float(gamma0) * v0v_ratio**(-float(q))
    return gamma_p

def e_thermal_integral(up_limit):
    npnts = 1000
    e_int = 0.0
    x = 0.0

    while x <= npnts:
        x += 1
        e_int += up_limit / npnts * (x * up_limit / npnts)**3.0 / (exp(x * up_limit / npnts) - 1)
    return e_int

def p_thermal(debye_t0, gamma0, q, n, t, v0v_ratio, v0):
    """Calculate P_thermal"""
    gamma_pp = gamma_p(gamma0, q, v0v_ratio)
    debye_tt = debye_t(debye_t0, gamma0, gamma_pp, q)
    et = 9 * n * R * t / (debye_tt / t)**3.0 * e_thermal_integral(debye_tt / t)
    e0 = 9 * n * R * 300.0 / (debye_tt / 300.0)**3.0 * e_thermal_integral(debye_tt / 300.0)

    p_thermal = gamma_pp / (v0 / v0v_ratio) * (et - e0)
    return p_thermal

def p_mgd(b0, b0p, v0, gamma0, q, theta0, n, pressure, temperature):
    """Calculate Pressure from Mie-Gruneisen-Deye EOS"""
    v0vratio = v0 / v_bm(b0, b0p, v0, pressure)
    p_mgd = pressure + p_thermal(theta0, gamma0, q, n, temperature, v0vratio, v0)
    return p_mgd

def v_mgd(b0, b0p, v0, gamma0, q, theta0, n, pressure, temperature):
    """calculate volume from Mie-Gruneisen-Debye EOS"""
    # iteration settings
    tol = 1e-9
    nmax = 1e4
    i = 1

    v0v_ratio = (1 + pressure * b0p / b0)**(1 / b0p)
    pressure_guess = bm_eos(b0, b0p, v0v_ratio)
    # iteration until p == p_init
    while n < nmax and abs(pressure_guess - pressure) > tol:
        pressure_guess = bm_eos(b0, b0p, v0v_ratio) + \
          p_thermal(theta0, gamma0, q, n, temperature, v0v_ratio, v0)
        v0v_ratio -=  (pressure_guess - pressure) / (b0 + b0p*pressure) / v0v_ratio
        i += 1
        if i == nmax: v0v_ratio = 1 
    volume = v0 / v0v_ratio
    return volume