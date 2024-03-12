#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 01:04:49 2024

@author: yan
"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import CoolProp as CP
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')

m1 = pd.read_csv("cfd_m1.csv", ",", skiprows=0)
model1 = pd.read_csv("result1.csv", ",", skiprows=0)


# m1 = pd.concat([m1,  ])
# axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


# fig 1
fig1 = plt.figure( dpi=300)
lw = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(model1.iloc[:,2] ,model1.iloc[:,-1], 'k+', lw=lw, label="MODEL")
axes.plot(m1.iloc[:,0] ,m1.iloc[:,1], 'ko', lw=lw, label="CFD")
# axes.plot(m1.iloc[:,0] ,m1.iloc[:,2], 'bo', lw=lw, label="CFD $P_2/P_t$")

# axes.set_xlim([0.5, 2])
# axes.set_ylim([0,0.7])
axes.set_xlabel('$X/D_e$',fontsize=12)
axes.set_ylabel('$P_a/P_t$',fontsize=12) 
# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("vv_nn_position.pdf")

# fig 2
fig2 = plt.figure( dpi=300)
lw = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(m1.iloc[:,1] ,m1.iloc[:,2]/m1.iloc[:,1], 'ko', lw=lw)

# axes.set_xlim([0.5, 2])
# axes.set_ylim([0,0.7])
axes.set_xlabel('$P_a/P_t$',fontsize=12)
axes.set_ylabel('$P_2/P_a$',fontsize=12) 
# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
# axes.legend(loc=0) # 

fig2.savefig("vv_nn_diff.pdf")
