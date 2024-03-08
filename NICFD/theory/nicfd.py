#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:56:12 2024

Find P2/Pt vs X/De for NICFD from input Pt.Tt, De/D or Me

@author: yan
"""

import CoolProp as CP
import math
import numpy as np
import pandas as pd
from RK4 import rk4

"""
0. fluid property
"""
fluidname = "nitrogen"
print("Fluid name:", fluidname)
R = CP.CoolProp.PropsSI("gas_constant",fluidname)
print("universal gas constant:  J/mol/K", R)
W = CP.CoolProp.PropsSI("molar_mass",fluidname)
print("molar mass: kg/mol", W)
Rs = R/W
print("spefific ags constant: J/Kg/K", Rs)
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)
print("critical temperature[K]:", Tc)
Pc =  CP.CoolProp.PropsSI("pcrit",fluidname)
print("critical pressure[Pa]:", Pc)

"""
1. input total conditions
"""

Pt = 8e5# total pressure
Tt = 337
s1 = CP.CoolProp.PropsSI('Smass','P',Pt,'T',Tt,fluidname) 
ht1 = CP.CoolProp.PropsSI('Hmass','P',Pt,'T',Tt,fluidname) 

"""
2. compute isentropic relationship
"""

n1 = 100
p = np.linspace(Pt*0.1,Pt,n1) # 
p = pd.Series(p)
h = np.zeros(p.size) # enthalpy
u = np.zeros(p.size) # velocity
c = np.zeros(p.size) # sound speed
m = np.zeros(p.size) # Mach number
t = np.zeros(p.size) # temperature


for i in p.index:
    t[i] = CP.CoolProp.PropsSI('T','P',p[i],'Smass',s1, fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','T',t[i],'P',p[i],fluidname)
    u[i] = math.sqrt(abs(2*(ht1-h[i])))
    c[i] = CP.CoolProp.PropsSI('A','P',p[i],'T',t[i],fluidname) 
    m[i] = u[i]/c[i] 

"""
3. find sonic condition, assume A* = 1
"""
print("index for sonic condition:",np.argmin(abs(m-1)))
Ms = m[np.argmin(abs(m-1))]
Ts = t[np.argmin(abs(m-1))]
Ps = p[np.argmin(abs(m-1))]
print("sonic P,T,M: ", Ps, Ts, Ms)
ds = CP.CoolProp.PropsSI('Dmass','T',Ts,'P',Ps,fluidname)
cs = c[np.argmin(abs(m-1))]
us = u[np.argmin(abs(m-1))]

"""
4. For Dr = De/Ds, find exit condition using constant mass flow
"""
De = 6.5
Ds = 6
Dr = De/Ds
n2 = 100
p = np.linspace(Ps*0.1,Ps,n2) # 
p = pd.Series(p)
h = np.zeros(p.size) 
d = np.zeros(p.size) 
u = np.zeros(p.size) # velocity
diff = np.zeros(p.size)
for i in p.index:
    d[i] = CP.CoolProp.PropsSI('Dmass','P',p[i],'Smass',s1, fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','Dmass',d[i],'P',p[i],fluidname)
    u[i] = math.sqrt(abs(2*(ht1-h[i])))
    diff[i] = ds*us - d[i]*u[i]*(Dr**2)
    
print("index for exit condition:",np.argmin(abs(diff)))
de = d[np.argmin(abs(m-1))]
Pe = p[np.argmin(abs(m-1))]
Te = CP.CoolProp.PropsSI('T','Dmass',de,'P',Pe,fluidname)
ue = u[np.argmin(abs(m-1))]
ce = CP.CoolProp.PropsSI('A','P|gas',Pe,'T',Te,fluidname)
Me = ue/ce 
print("exit P,T,M: ", Pe, Te, Me)    

"""
5.Find X/De by solving MOC
"""

dx = De/100 # Delta x
n3 = 200 
k = np.zeros(n3-1) # number of inetraion, n-1
k = pd.Series(k)
# initialization
x = np.zeros(n3)
x = pd.Series(x)
x[0] = 0 # 
y = np.zeros(n3)
y[0] = -De/2 #  
dy = np.zeros(n3)
dy[0] = 0 #  
M = np.zeros(n3)
M[0] = Me # Mach number

vvv,ttt,mmm,nununu = rk4(1/ds, 10/ds, Ts, 1, 0, 200)

mu = np.zeros(n3) 
mu[0] = np.arcsin(1/Me) # Mach angle
nu = np.zeros(n3) 
nu[0] = nununu[np.argmin(abs(mmm-Me))] # Prandtl Meyer angle
print("index for nu0:",np.argmin(abs(mmm-Me)))
print("initial nu:", nu[0])
theta = np.zeros(n3)
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
    M[i+1] = mmm[np.argmin(abs(nununu-nu[i+1]))]
    mu[i+1] = np.arcsin(1/M[i+1])
    if y[i+1]>0:
        M1 = M[i+1]
        xx = x[i]/De
        p = np.linspace(Pe*0.05,Pe,n2) # 
        p = pd.Series(p)
        h = np.zeros(p.size) 
        d = np.zeros(p.size) 
        u = np.zeros(p.size) # velocity
        diff = np.zeros(p.size)
        for i in p.index:
            c[i] = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s1,fluidname)
            h[i] = CP.CoolProp.PropsSI('Hmass','P',p[i],'Smass',s1,fluidname)
            u[i] = M1*c[i]
            diff[i] = ht1 - h[i] - 0.5*u[i]*u[i] 
        print("index for P1:" , np.argmin(abs(diff)) )
        P1 = p[np.argmin(abs(diff))]
        T1 = CP.CoolProp.PropsSI('T','P',P1,'Smass',s1,fluidname)
        print("centerline condition; X/De, M, P[Pa], T[k] " , xx, M1, P1, T1 )
        break

"""
6.Find post-shock states
"""