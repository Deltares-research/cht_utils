# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:26:53 2022

@author: Roel de Goede
"""

from .disper import disper
import numpy as np
from scipy import interpolate as intp
    
def deshoal(hm0,Tp,d_profile,d_BC):
    cg_profile = wavecelerity(Tp,d_profile)
    cg_BC = wavecelerity(Tp,d_BC)
    
    hm0_deshoal = hm0 * np.sqrt(cg_profile/cg_BC)
    return hm0_deshoal
    
def wavecelerity(Tp, d):
    g= 9.81
    k   = disper(2*np.pi/Tp, d, g)
    n   = 0.5*(1+2*k*d/np.sinh(2*k*d))
    c   = g*Tp/(2*np.pi)*np.tanh(k*d)
    cg  = n*c
    return cg