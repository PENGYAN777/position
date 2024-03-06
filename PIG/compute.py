#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 01:10:23 2024

Compute P/Pe for PIG, in terms of sonic and supersonic nozzle

@author: yan
"""

import math
import numpy as np
import pandas as pd
from PM_function import NufromM, MformNu

"""
Sonic model
"""

# Total condition, taken from experiment of under-expanded jet with nitrogen, t=5s.
Pt = 4.57
Tt = 311.35
gamma = 1.4 # cp/cv

# sonic condition
mue = 85
Me = 1/math.sin(math.radians(mue))
Pe = Pt*(1+(gamma-1)/2*Me*Me)**(-gamma/(gamma-1))
Te = Tt*(1+(gamma-1)/2*Me*Me)**(-1)


# find P(x=0) along centerline using MOC
De = 1 # diameter of nozzle exit
dx = De/100 # Delta x 
n = 200 
k = np.zeros(n-1) # number of inetraion, n-1
k = pd.Series(k)

################## find P,M at y=0 for different theta
Theta = np.linspace(0,1,5)*math.pi/180
Theta = pd.Series(Theta)

xx = np.zeros(Theta.size) # x with y=0
MM = np.zeros(Theta.size) # M with y=0
PP = np.zeros(Theta.size) # P/Pe with y=0


for j in Theta.index:
    print("j= ", j)
    # initilization
    x = np.zeros(n)
    x = pd.Series(x)
    x[0] = 0 # 
    y = np.zeros(n)
    y[0] = De/2 #  
    M = np.zeros(n)
    M[0] = Me # Mach number
    mu = np.zeros(n) 
    mu[0] = mue # list of Mach angle mu
    nu = np.zeros(n) 
    nu[0] = NufromM(M[0]) # list Prandtl Meyer function
    theta = np.zeros(n)
    theta[0] = Theta[j] # theta
    for i in k.index:
        dy = dx*math.tan(theta[i]-mu[i])
        x[i+1] = x[i] + dx
        y[i+1] = y[i] + dy
        theta[i+1] = np.arctan(dy/dx)
        nu[i+1] = nu[i] -theta[i+1] + theta[i] 
        M[i+1] = MformNu(nu[i+1])
        mu[i+1] = np.arcsin(1/M[i+1])
    x[:] = x[::-1]     
    y[:] = y[::-1] # reverse y so that once find y>0, stop 
    M[:] = M[::-1]  
    for i in x.index:
        if y[i]>0:
            MM[j] = M[i]
            PP[j] = Pt*(1+(gamma-1)/2*MM[j]*MM[j])**(-gamma/(gamma-1))/Pe
            xx[j] = x[i]/De
            break
    # MM[j] = M[np.argmin(abs(y))] # the Mach number at centerline
    # PP[j] = Pt*(1+(gamma-1)/2*MM[j]*MM[j])**(-gamma/(gamma-1))/Pe
    # xx[j] = x[np.argmin(abs(y))]/De
    