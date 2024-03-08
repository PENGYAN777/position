#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:16:39 2024

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

sonic = pd.read_csv("sonic.csv", ",", skiprows=0)
m = pd.read_csv("paper_m.csv", ",", skiprows=0)

# fig 1
fig1 = plt.figure( dpi=300)
lwh = 2
axes = fig1.add_axes([0.15, 0.15, 0.7, 0.7]) #size of figure
axes.plot(sonic.iloc[:,2] , sonic.iloc[:,3], 'k+', lw=lwh, label="Theoretical solution")
axes.plot(m.iloc[:,0] , m.iloc[:,1], 'ko', lw=lwh, label="PhilipL 1948")

# compute diff
ax2 = axes.twinx()
idx = np.zeros(m.iloc[:,0].size) # 
diff = np.zeros(m.iloc[:,0].size) # diff between cfd and ex
for i in m.iloc[:,0].index:
    idx[i] = np.argmin(abs(m.iloc[i,0]-sonic.iloc[:,2])) # 
    diff[i] = (m.iloc[i,1] - sonic.iloc[int(idx[i]),3] )/(m.iloc[i,1])*100
    
ax2.plot(m.iloc[:,0], abs(diff) , 'b*', lw=lwh)    
ax2.set_ylim([0,20])
# ax2.set_ylabel('$\Delta_{P/P_t}$(%)',fontsize=12) 

axes.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axes.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

axes.set_xlim([0,2])
axes.set_ylim([1,4])
axes.set_xlabel('$X/D_e$',fontsize=12)
axes.set_ylabel('$M$',fontsize=12) 
# axes.set_title('$P/P_c$ along y$=0.3 mm$',fontsize=14)
axes.legend(loc=0) # 

fig1.savefig("sonic_vv_m.pdf")