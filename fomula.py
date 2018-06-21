# -*- coding: utf-8 -*-
"""
Created on Tue Feb 06 04:35:03 2018

@author: Administrator
"""
from __future__ import absolute_import, division, print_function
import numpy as np
from const import (dict, mu_func, si_func)

def RK(muvar, N, T0, h0):    
    # The function of Runge-Kutta numerical schemes
    T = np.zeros(N)
    h = np.zeros(N)
    t = np.zeros(N)    
    T[0] = T0/7.5 # non-dimensionalised
    h[0] = h0/150 # non-dimensionalised
    t[0] = 0.
    
    for i in range(0,N-1):
        t[i+1] = t[0] + dict['dt']*(i+1)
        miu = mu_func(muvar,t[i]) # The coupling coefficient, miu, dimensionless
        pesi = si_func(t[i]) # The parameters for calculating the additional wind stress forcing
        k1, l1 = f(T[i],h[i],miu,pesi)
        k2, l2 = f(T[i] + dict['dt']*k1/2, h[i] + dict['dt']*l1/2,miu,pesi)
        k3, l3 = f(T[i] + dict['dt']*k2/2, h[i] + dict['dt']*l2/2,miu,pesi)
        k4, l4 = f(T[i] + dict['dt']*k3, h[i] + dict['dt']*l3,miu,pesi)        
        T[i+1] = T[i] + (k1 + 2*k2 + 2*k3 + k4)*dict['dt']/6.
        h[i+1] = h[i] + (l1 + 2*l2 + 2*l3 + l4)*dict['dt']/6.
    return T, h, t
    
def f(T0, h0, miu, pesi):
    # The function of Fei-Fei Jin's ENSO recharge model
    R = dict['gama']*dict['b0']*miu - dict['c']
    # R collectively describes the Bjerknes positive feedback process
    T = R*T0 + dict['gama']*h0 - dict['epsilon']*(h0 + dict['b0']*miu*T0)**3 \
        + dict['gama']*pesi
    h = -dict['r']*h0 - dict['alpha']*miu*dict['b0']*T0 - dict['alpha']*pesi
    
    return T, h
    
