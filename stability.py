import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import cmath
from const import (dict, mu_func, si_func)

def ampli(mu):
    # The function for calculating the amplification factors
    b=dict['b0']*mu
    # b relates a stronger thermocline gradient to stronger easterly wind stress
    R=dict['gama']*b-dict['c']
    # R collectively describes the Bjerknes positive feedback process
    d=R
    e=dict['gama']
    f=-dict['alpha']*b
    g=-dict['r']

    A=np.empty([2, 2], float) # define the coefficient matrix
    ampl = np.zeros(2)
    v = np.zeros(2)
    A=np.array([[d,e],[f,g]])
    eigenvalues,eigenvectors = np.linalg.eig(A)

    v[0] = pl.real(dict['dt']*eigenvalues[0])
    v[1] = pl.real(dict['dt']*eigenvalues[1])

    ampl[0] = 1+v[0]+0.5*v[0]**2+1./6.*v[0]**3+1./24.*v[0]**4
    ampl[1] = 1+v[1]+0.5*v[1]**2+1./6.*v[1]**3+1./24.*v[1]**4
    # calculate the amplification factors of 4th-order Runge-Kutta scheme

    return ampl[0],ampl[1]


