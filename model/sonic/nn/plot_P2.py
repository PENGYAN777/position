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

result1 = pd.read_csv("result1.csv", ",", skiprows=0)
result = pd.concat([result1,  ])
# axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
gamma = 1.4

# fig 1
fig1 = plt.figure( dpi=300)
lw = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(result.iloc[:,2] ,result.iloc[:,-1], 'ko', lw=lw, label="$P_1/P_t$")
axes.plot(result.iloc[:,2] ,result.iloc[:,-2], 'bo', lw=lw, label="$P_2/P_t$")

# axes.set_xlim([0.5, 2])
# axes.set_ylim([0,0.7])
axes.set_xlabel('$X/D_e$',fontsize=12)
axes.set_ylabel('$P/P_t$',fontsize=12) 
# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("example_nn_p.pdf")

# fig 1
fig2 = plt.figure( dpi=300)
lw = 2
axes = fig2.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(result.iloc[:,2] ,result.iloc[:,3], 'ko', lw=lw, label="$M_1$")
axes.plot(result.iloc[:,2] ,result.iloc[:,4], 'bo', lw=lw, label="$M_2$")

# axes.set_xlim([0.5, 2])
# axes.set_ylim([0,0.7])
axes.set_xlabel('$X/D_e$',fontsize=12)
axes.set_ylabel('$M$',fontsize=12) 
# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
axes.legend(loc=0) # 

fig2.savefig("example_nn_mach.pdf")