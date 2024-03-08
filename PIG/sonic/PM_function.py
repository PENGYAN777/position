#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 03:41:21 2024

@author: user
"""

import math
import numpy as np
import pandas as pd

def NufromM(M):
    """
    Input Mach number and output nu
    """
    g = 1.4 # cp/cv
    f1 = math.sqrt((g+1)/(g-1))
    f2 = math.sqrt((g-1)/(g+1)*(M*M-1))
    f3 = math.sqrt(M*M-1)
    nu = f1*np.arctan(f2) - np.arctan(f3)
    return nu


def MformNu(nu):
    n = 500
    M = np.linspace(1,10,n)
    M = pd.Series(M)
    diff = np.zeros(M.size)
    for i in M.index:
        diff[i] = NufromM(M[i]) - nu
    m = M[np.argmin(abs(diff))]
    return m