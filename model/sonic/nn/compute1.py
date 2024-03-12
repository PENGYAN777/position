#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:56:12 2024

Find P2/Pt vs X/De for NICFD from input Pt.Tt, De/D or Me
Use d for denisty and D for diameter

@author: yan
"""

import CoolProp as CP
import math
import numpy as np
import pandas as pd
from function import *
import time

# Start timer
start_time = time.perf_counter()

"""
0. fluid property
"""

g = 1.4 # gamma = cp/cv


"""
1. input total conditions
"""
# non-dimensional
Pt = 1
Tt = 1

"""
2. find sonic condition, assume A* = 1.001
"""

Ms = 1.001
Ps = PPtFromM(Ms,g)*Pt

"""
3. Loop for different Mach 
"""
# prepare NufromM
Mi = np.linspace(Ms, 10, 100)
Mi = pd.Series(Mi)
Nui = np.zeros(Mi.size) 
for i in Mi.index:
    Nui[i] = NuFromM(Mi[i],g)

# range of Mach
Me = np.linspace(Ms, 5, 20)
Me = pd.Series(Me)
Pe = np.zeros(Me.size)  
De = 1 # Diameter of nozzle exit
xx = np.zeros(Me.size)  
M1 = np.zeros(Me.size)  
M2 = np.zeros(Me.size)  
P1 = np.zeros(Me.size)  
P2 = np.zeros(Me.size)  

dx = De*0.005 # Delta x
n1 = 600  # large n for large M
print("max X/De: ", n1*dx)
k = np.zeros(n1-1) # number of inetraion, n-1
k = pd.Series(k)
for j in Me.index:
    print("j = ", j)
    """
    3.1 Find Pe, Te for different Me
    """
    Pe[j] = PPtFromM(Me[j],g)*Pt
    """
    3.2 Find X/De by solving MOC
    """
    # initialization
    x = np.zeros(n1)
    x = pd.Series(x)
    x[0] = 0 # 
    y = np.zeros(n1)
    y[0] = -De/2 #  
    dy = np.zeros(n1)
    dy[0] = 0 #  
    M = np.zeros(n1)
    M[0] = Me[j] # Mach number
    mu = np.zeros(n1) 
    mu[0] = np.arcsin(1/M[0]) # Mach angle
    nu = np.zeros(n1) 
    nu[0] = NuFromM(M[0], g) # Prandtl Meyer angle
    theta = np.zeros(n1)
    theta[0] = 0 # flow angle
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
        M[i+1] = Mi[np.argmin(abs(Nui-nu[i+1]))]
        mu[i+1] = np.arcsin(1/M[i+1])
        if y[i+1]>0:
            M1[j] = M[i+1]
            xx[j] = x[i]/De
            P1[j] = PPtFromM(M1[j],g)*Pt
            break
    """
    3.3.Find post-shock states
    """
    M2[j] = M2FromM1(M1[j],g)
    P2[j] = P1[j]*P2P1FromP1(M1[j],g)
    print("------------------------------------------------------------")
    print("post-shock condition; M1, M2,P1/Pt, P2/Pt" ,M1[j] , M2[j],P1[j]/Pt, P2[j]/Pt )
    print("------------------------------------------------------------")


"""
4. write results into csv file
"""
pd.DataFrame(xx).to_csv('result1.csv', index_label = "Index", header  = ['X/De'])
data = pd.read_csv("result1.csv", ",")
D =pd.DataFrame({'M1': M1, 'M2': M2, 'P1/Pt': P1/Pt, 'P2/Pt': P2/Pt,})
newData = pd.concat([data, D], join = 'outer', axis = 1)
newData.to_csv("result1.csv")










"""
xx. output computational time
"""
# End timer
end_time = time.perf_counter()


# Calculate elapsed time
elapsed_time = end_time - start_time
print("Elapsed time: ", elapsed_time)