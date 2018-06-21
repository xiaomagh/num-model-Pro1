# -*- coding: utf-8 -*-
"""
student number:25806676
"""
from __future__ import division
import random
import numpy as np

dict = {
# High-end value of the coupling parameter, b0, dimensionless
   'b0': 2.5,
# The constant specifies the feedback of the thermocline gradient 
# on the SST gradient, gama, dimensionless
   'gama': 0.75,
# The damping rate of SST anomalies, c, dimensionless
   'c': 1.,
# The parameter represents the damping of the upper ocean heat content, r, dimensionless
   'r': 0.25,
# The parameter relates enhanced easterly wind stress to the 
#recharge of ocean heat content, alpha, dimensionless
   'alpha': 0.125,
# The parameter varies the degree of nonlinearity, epsilon, dimensionless
   'epsilon': 0.,
# Time step after non-dimensionlised
   'dt': 0.5,
# The coupling coefficient, miu, dimensionless
   'miu': 2./3.,
# The parameters for calculating the additional wind stress forcing
   'fann': 0.02, 'fran':0.2,
# The parameters for calculating the coupling coefficient, miu
   'miu0': 0.75, 'miua': 0.2,
   }

type(dict)

# The function of wind stress forcing of the system, pesi
def si_func(t):
    pesi = dict['fann']*np.cos(2*np.pi*t/6.) + dict['fran']*random.uniform(-1,+1)
    return pesi 
# The function of coupling coefficient, miu
def mu_func(muvar,t):
    if(muvar):
        m = dict['miu0']*(1 + dict['miua']*np.cos(2*np.pi*t/6. - 5*np.pi/6.))
    else:
        m = dict['miu']
    return m

            
