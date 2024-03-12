#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 17:10:25 2023

@author: yan
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import CoolProp as CP



"""
Isentropic relation
"""

def PPtFromM(M,g):
    """

    Parameters
    ----------
    M : Mach number
    g : cp/cv

    Returns P/Pt
    -------
    """
    f1 = 1+(g-1)/2*M*M
    f2 = -g/(g-1)
    PPt = f1**(f2)
    return PPt

def NuFromM(M,g):
    """

    Parameters
    ----------
    M : Mach number
    g : cp/cv

    Returns Prandtl Meyer angle
    -------

    """
    f1 = (g+1)/(g-1)
    f2 = (g-1)/(g+1)*(M*M-1)
    f3 = M*M-1
    nu = math.sqrt(f1)*math.atan(math.sqrt(f2)) - math.atan(math.sqrt(f3))
    return nu

"""
Normal shock relation
"""

def M2FromM1(M1,g):
    """

    Parameters
    ----------
    M : Mach number
    g : cp/cv

    Returns M2
    -------
    """
    M2 = ((g-1)*M1*M1+2)/(2*g*M1*M1-g+1)
    M2 = math.sqrt(M2)
    return M2
    
def P2P1FromP1(M1,g):
    """

    Parameters
    ----------
    M : Mach number
    g : cp/cv

    Returns P2/P1
    -------
    """
    P2 = 2*g*M1*M1/(g+1) - (g-1)/(g+1)
    return P2
    