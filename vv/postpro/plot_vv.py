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

cfd_nn = pd.read_csv("cfd_nn.csv", ",", skiprows=0)
model_nn = pd.read_csv("model_nn.csv", ",", skiprows=0)

# axes.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


# fig 1
fig1 = plt.figure( dpi=300)
lw = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(model_nn.iloc[:,2] ,model_nn.iloc[:,-3], 'ko', lw=lw, label="$Model$")
axes.plot(cfd_nn.iloc[:,1] , 1/cfd_nn.iloc[:,0], 'bo', lw=lw, label="$CFD$")

# axes.set_xlim([0.4, 1.5])
# axes.set_ylim([0,1.0])
axes.set_xlabel('$X/D_e$',fontsize=12)
axes.set_ylabel('$P_a/P_t$',fontsize=12) 
# axes.set_title('$P/P_t$ along nozzle centerline',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("vv_nn.pdf")