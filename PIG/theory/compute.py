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
Ms = 1/math.sin(85*math.pi/180 )
Me = np.linspace(1,4,20)
Me = pd.Series(Me)
mue = np.arcsin(1/Me)
# Te = Tt*(1+(gamma-1)/2*Me*Me)**(-1)


# find P(x=0) along centerline using MOC
De = 1 # diameter of nozzle exit
dx = De/100 # Delta x 
n = 200 
k = np.zeros(n-1) # number of inetraion, n-1
k = pd.Series(k)

################## find P,M at y=0 for different theta
# Theta = np.linspace(0,1,1)*math.pi/180 # degree to radius
# Theta = pd.Series(Theta)

xx = np.zeros(mue.size) # x with y=0
MM = np.zeros(mue.size) # M with y=0
PP = np.zeros(mue.size) # P/Pe with y=0



for j in mue.index:
    # initilization
    print("j= ", j, "Me= ", Me[j])
    Pe = Pt*(1+(gamma-1)/2*Ms*Ms)**(-gamma/(gamma-1))
    x = np.zeros(n)
    x = pd.Series(x)
    x[0] = 0 # 
    y = np.zeros(n)
    y[0] = -De/2 #  
    dy = np.zeros(n)
    dy[0] = 0 #  
    M = np.zeros(n)
    M[0] = Me[j] # Mach number
    mu = np.zeros(n) 
    mu[0] = mue[j] # list of Mach angle mu
    nu = np.zeros(n) 
    nu[0] = NufromM(M[0]) # list Prandtl Meyer function
    theta = np.zeros(n)
    theta[0] = 0 # theta
    for i in k.index:
        dy[i] = dx*math.tan(theta[i]+mu[i])
        x[i+1] = x[i] + dx
        y[i+1] = y[i] + dy[i]
        dtheta = np.linspace(0,1,100)*math.pi/180 # degree to radius
        dtheta = pd.Series(dtheta)
        diff = np.zeros(dtheta.size)
        for ii in dtheta.index:
            t1 = theta[i]
            t2 = theta[i] + dtheta[ii] # t1,t2 temporary variable
            diff[ii] = math.cos(t2) - math.cos(t1) - dy[i]/dx*(math.sin(t2) - math.sin(t1)) 
        theta[i+1] = theta[i] + dtheta[np.argmin(abs(diff))]
        nu[i+1] = nu[i]  + theta[i+1] - theta[i] 
        M[i+1] = MformNu(nu[i+1])
        mu[i+1] = np.arcsin(1/M[i+1])
        if y[i+1]>0:
            MM[j] = M[i+1]
            PP[j] = Pt*(1+(gamma-1)/2*MM[j]*MM[j])**(-gamma/(gamma-1))/Pe
            xx[j] = x[i]/De
            break

# write results
pd.DataFrame(xx).to_csv('sonic.csv', index_label = "Index", header  = ['X/De']) 
data = pd.read_csv("sonic.csv", ",")
D =pd.DataFrame({'M': MM, 'P/Pe': PP})
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("sonic.csv")
