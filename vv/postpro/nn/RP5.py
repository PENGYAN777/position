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
from RK4 import rk4
import time

# Start timer
start_time = time.perf_counter()

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

# interasted variables need to write into csv file
xx = []
P2 = []
T2 = []
M2 = []
P1 = []

"""
1. input total conditions
"""

Pt = 5e5# total pressure
Tt = 337
s1 = CP.CoolProp.PropsSI('Smass','P',Pt,'T',Tt,fluidname) 
ht1 = CP.CoolProp.PropsSI('Hmass','P',Pt,'T',Tt,fluidname) 
Pa = 1e5 # ambient pressure

"""
2. compute isentropic relationship
"""

n1 = 100
p = np.linspace(Pt,Pt*0.1,n1) # 
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
    if abs(m[i]-1)<0.01:
        break
"""
3. find sonic condition, assume A* = 1
"""
print("index for sonic condition:",np.argmin(abs(m-1)))
Ms = m[i]
Ts = t[i]
Ps = p[i]
print("sonic P,T,M: ", Ps, Ts, Ms)
ds = CP.CoolProp.PropsSI('Dmass','T',Ts,'P',Ps,fluidname)
cs = c[i]
us = u[i]

# use sonic condition to compute Prandtl Meyer angle
vvv,ttt,mmm,nununu = rk4(1/ds, 30/ds, Ts, 1, 0, 200)

"""
4. For Dr = De/Ds, find nozzle exit condition using constant mass flow
"""
De = 6.5
Ds = 6
Dr = De/Ds
n2 = 100
p = np.linspace(Ps,Ps*0.1,n2) # 
p = pd.Series(p)
h = np.zeros(p.size) 
d = np.zeros(p.size) 
u = np.zeros(p.size) # velocity
diff = np.zeros(p.size)
for i in p.index:
    d[i] = CP.CoolProp.PropsSI('Dmass','P',p[i],'Smass',s1, fluidname) 
    h[i] = CP.CoolProp.PropsSI('Hmass','Dmass',d[i],'P',p[i],fluidname)
    u[i] = math.sqrt(abs(2*(ht1-h[i])))
    diff[i] = (ds*us - d[i]*u[i]*(Dr**2))/ds/us
    if abs(diff[i])<0.01:
        break
# print("index for exit condition:",i)
de = d[i]
pe = p[i]
Te = CP.CoolProp.PropsSI('T','Dmass',de,'P',pe,fluidname)
ue = u[i]
ce = CP.CoolProp.PropsSI('A','P|gas',pe,'T',Te,fluidname)
me = ue/ce 
# print("exit P,T,M: ", Pe, Te, Me)    
print("Emperical: ", math.sqrt(1.4/2*pe/Pa*me*me))
"""
5. Loop for different Me
"""
Me = np.linspace(me,4, 10)
Me = pd.Series(Me)
print("max Me < max mmm !!!  : ", Me.iloc[-1], mmm[-1])
Pe = np.zeros(Me.size) 
for j in Me.index:
    print("j = ", j)
    """
    5.0 Find Pe, Te for different Me
    """
    p = np.linspace(pe,pe*0.01,500) # 
    p = pd.Series(p)
    h = np.zeros(p.size) 
    c = np.zeros(p.size) 
    u = np.zeros(p.size) 
    diff = np.zeros(p.size)
    for i in p.index:
        c[i] = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s1,fluidname)
        u[i] = Me[j]*c[i]
        h[i] = CP.CoolProp.PropsSI('Hmass','P',p[i],'Smass',s1,fluidname)
        diff[i] = (ht1 - h[i] - 0.5*u[i]*u[i])/ht1 
        if abs(diff[i])<0.01:
            Pe[j] = p[i]
            break    
    """
    5.1 Find X/De by solving MOC
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
    M[0] = Me[j] # Mach number
    
    mu = np.zeros(n3) 
    mu[0] = np.arcsin(1/M[0]) # Mach angle
    nu = np.zeros(n3) 
    nu[0] = nununu[np.argmin(abs(mmm-M[0]))] # Prandtl Meyer angle
    # print("index for nu0:",np.argmin(abs(mmm-Me)))
    # print("initial nu:", nu[0])
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
            xx.append(x[i]/De)
            p = np.linspace(Pe[j],Pe[j]*0.1,n2) # 
            p = pd.Series(p)
            h = np.zeros(p.size) 
            c = np.zeros(p.size) 
            u = np.zeros(p.size) 
            diff = np.zeros(p.size)
            for i in p.index:
                c[i] = CP.CoolProp.PropsSI('A','P',p[i],'Smass',s1,fluidname)
                h[i] = CP.CoolProp.PropsSI('Hmass','P',p[i],'Smass',s1,fluidname)
                u[i] = M1*c[i]
                diff[i] = (ht1 - h[i] - 0.5*u[i]*u[i])/ht1 
                if abs(diff[i])<0.01:
                    break
            # print("index for P1:" , i )
            P1.append(p[i])
            T1 = CP.CoolProp.PropsSI('T','P',p[i],'Smass',s1,fluidname)
            d1 = CP.CoolProp.PropsSI('Dmass','P', p[i],'T',T1,fluidname)
            u1 = u[i]
            # print("centerline condition; X/De, M, P[Pa], T[k] " , xx, M1, P1, T1 )
            break
    
    """
    5.2.Find post-shock states
    """
    n4 = 1000
    p = np.linspace(P1[j]*1.1,P1[j]*(6*Me[j] - 4),n4) # must choose reasonable range
    p = pd.Series(p)
    h = np.zeros(p.size) 
    d = np.zeros(p.size) 
    u = np.zeros(p.size) 
    diff = np.zeros(p.size)
    for i in p.index:
        u[i] = (P1[j]+d1*u1*u1-p[i])/d1/u1
        if u[i]>0:
            d[i] =  d1*u1/u[i]
            h[i] = CP.CoolProp.PropsSI('Hmass','P',p[i],'Dmass',d[i],fluidname) 
            diff[i] = (ht1 - 0.5*u[i]*u[i] - h[i])/ht1
            if abs(diff[i])<0.01:
                break
    P2.append(p[i])
    d2 = d[i]
    T2.append( CP.CoolProp.PropsSI('T','P',p[i],'Dmass',d2,fluidname)  )
    c2 = CP.CoolProp.PropsSI('A','P',p[i],'Dmass',d2,fluidname) 
    u2 = u[i]
    M2.append( u2/c2 ) 
    print("------------------------------------------------------------")
    print("post-shock condition; index, M1, M2, P2/Pt" ,  i,Me[j] ,u2/c2, p[i]/Pt )
    print("------------------------------------------------------------")
    if p[i] < Pa:
        print("------------------------------------------------------------")
        print("Find P2=Pa; index, Pa, P2" ,  j, Pa, p[i] )
        print("------------------------------------------------------------")
        break

"""
6. write results into csv file
"""
# convert list to array
# xx = np.array(xx)
# P2 = np.array(P2)
# P1 = np.array(P1)
# T2 = np.array(T2)
# M2 = np.array(M2)

# pd.DataFrame(xx).to_csv('result.csv', index_label = "Index", header  = ['X/De'])
# data = pd.read_csv("result.csv", ",")
# D =pd.DataFrame({'Me': Me, 'M2': M2, 'P2/Pt': P2/Pt, 'T2/Tt': T2/Tt,'P1/Pt': P1/Pt,})
# newData = pd.concat([data, D], join = 'outer', axis = 1)
# newData.to_csv("result.csv")



"""
xx. output computational time
"""
# End timer
end_time = time.perf_counter()


# Calculate elapsed time
elapsed_time = end_time - start_time
print("Elapsed time: ", elapsed_time)