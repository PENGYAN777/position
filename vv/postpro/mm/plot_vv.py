#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 02:12:47 2024

@author: yan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os
import CoolProp as CP
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')

cfd = pd.read_csv("cfd.csv", ",", skiprows=0)
model = pd.read_csv("model.csv", ",", skiprows=0)

# axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


# fig 1
fig1 = plt.figure( dpi=300)
lw = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(model.iloc[:,0] ,model.iloc[:,1], 'ko', lw=lw, label="$Model$")
axes.plot(cfd.iloc[:,0] , cfd.iloc[:,1], 'bo', lw=lw, label="$CFD$")

axes.set_xlim([4, 10])
axes.set_ylim([0,2])
axes.set_xlabel('$P_t/P_a$',fontsize=12) 
axes.set_ylabel('$X/D_e$',fontsize=12)

# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("vv_mm.pdf")