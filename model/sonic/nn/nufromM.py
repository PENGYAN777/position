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

fluidname = "nitrogen"
Tc =  CP.CoolProp.PropsSI("Tcrit",fluidname)

# function to be solved
def f1(v,M):
    return  math.sqrt(M*M-1)/v/M/M

def f2(v,T,M):
    if abs(T-Tc)<0.01*Tc:
        T = 0.99*Tc
    g = CP.CoolProp.PropsSI('fundamental_derivative_of_gas_dynamics','T',T,'Dmass',1/v,fluidname ) 
    return  -M/v*(1 - g - 1/M/M)

def f3(v,T):
    if abs(T-Tc)<0.01*Tc:
        T = 0.99*Tc
    P = CP.CoolProp.PropsSI('P','T',T,'Dmass',1/v,fluidname ) 
    return  -T/CP.CoolProp.PropsSI('Cvmass','T',T,'Dmass',1/v,fluidname)* CP.CoolProp.PropsSI('d(P)/d(T)|Dmass','P',P,'T',T,fluidname )

def rk4(v0,vn,t0,m0,nu0,n):
    dv = (vn-v0)/n          # estep Delta v
    h = dv
    """
    initialization
    """
    v_arr = np.zeros(n + 1)   # create an array of zeros for spefific volume
    t_arr = np.zeros(n +1)    # create an array of zeros for Temperature
    m_arr = np.zeros(n + 1)   # create an array of zeros for Mach number
    nu_arr = np.zeros(n + 1)   # create an array of zeros for Prandtl Mayer function nu
    
    v_arr[0] = v0              # add initial value to array
    t_arr[0] = t0              # add initial value to array
    m_arr[0] = m0              # add initial value to array
    nu_arr[0] = nu0              # add initial value to array

    # Euler's method
    for i in range (1, n + 1):  
       v = v_arr[i-1]
       t = t_arr[i-1]
       m = m_arr[i-1]
       nu = nu_arr[i-1]
       # update nu
       k1 = h * (f1(v, m))
       k2 = h * (f1(v+h/2, m))
       k3 = h * (f1(v+h/2, m))
       k4 = h * (f1(v+h, m))
       k = (k1+2*k2+2*k3+k4)/6
       nu_arr[i] = nu + k
       # update M
       k1 = h * (f2(v, t,m))
       k2 = h * (f2(v+h/2,t+h/2,m+1/2*k1))
       k3 = h * (f2(v+h/2,t+h/2, m+1/2*k2))
       k4 = h * (f2(v+h,t+h, m+k3))
       k = (k1+2*k2+2*k3+k4)/6
       m_arr[i] = m + k
       # update T
       k1 = h * (f3(v, t))
       k2 = h * (f3(v+h/2,t+1/2*k1))
       k3 = h * (f3(v+h/2,t+1/2*k2))
       k4 = h * (f3(v+h,t+k3))
       k = (k1+2*k2+2*k3+k4)/6
       t_arr[i] = t + k
       # update v
       v_arr[i] = v + h
       
       # dnudv = math.sqrt(m*m-1)/v/m/m
       # dmdv = -m/v*(1 - g - 1/m/m)
       # dtdv = -t/CP.CoolProp.PropsSI('Cvmass','T',t,'Dmass',1/v,"Toluene")* CP.CoolProp.PropsSI('d(P)/d(T)|Dmass','P',p,'T',t,fluidname )
       # v_arr[i] = v + dv
       # t_arr[i] = t + dv*dtdv
       # m_arr[i] = m + dv*dmdv
       # nu_arr[i] = nu + dv*dnudv
       
    return v_arr,t_arr,m_arr,nu_arr
