#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 01:10:23 2024

Compute P/Pe for PIG, in terms of sonic and supersonic nozzle

@author: yan
"""

import math

"""
Sonic model
"""

# Total condition, taken from experiment of under-expanded jet with nitrogen, t=5s.
Pt = 4.57
Tt = 311.35
gamma = 1.4 # cp/cv

# sonic condition
mue = 85
Me = 1/math.sin(math.radians(mue))
Pe = Pt*(1+(gamma-1)/2*Me*Me)**(-gamma/(gamma-1))
Te = Tt*(1+(gamma-1)/2*Me*Me)**(-1)


# find P(x=0) along centerline using MOC
