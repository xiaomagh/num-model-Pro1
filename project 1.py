# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 19:29:53 2018

@author: Administrator
"""

def taska(r = 0.25, alpha = 0.125, b = 2.5, miu = 2./3., N = 121, dt = 0.5, \
          gama = 0.75, name_fig = 'TASK A'):

    # Set the variable lists
    T = np.zeros(N)
    h = np.zeros(N)
    t = np.zeros(N)
    # Define the Runge-Kutta scheme 
    k1 = np.zeros(N-1)
    k2 = np.zeros(N-1)
    k3 = np.zeros(N-1)
    k4 = np.zeros(N-1)
    l1 = np.zeros(N-1)
    l2 = np.zeros(N-1)
    l3 = np.zeros(N-1)
    l4 = np.zeros(N-1)

    T[0] = 1.125
    h[0] = 0.
    t[0] = 0.

    for i in range(1,N):
        t[i] = 0. + dt*i
    
        k1[i-1] = r*T[i-1] + gama*h[i-1]
        l1[i-1] = -r*h[i-1] - alpha*b*miu*T[i-1]
        k2[i-1] = r*(T[i-1] + (dt/2)*k1[i-1]) + gama*(h[i-1] + (dt/2)*l1[i-1])
        l2[i-1] = -r*(h[i-1] + (dt/2)*l1[i-1]) - \
                  alpha*b*miu*(T[i-1] + (dt/2)*k1[i-1])
        k3[i-1] = r*(T[i-1] + (dt/2)*k2[i-1]) + gama*(h[i-1] + (dt/2)*l2[i-1])
        l3[i-1] = -r*(h[i-1] + (dt/2)*l2[i-1]) - \
                  alpha*b*miu*(T[i-1] + (dt/2)*k2[i-1])
        k4[i-1] = r*(T[i-1] + dt*k3[i-1]) + gama*(h[i-1] + dt*l3[i-1])
        l4[i-1] = -r*(h[i-1] + dt*l3[i-1]) - alpha*b*miu*(T[i-1] + dt*k3[i-1])
    
        T[i] = T[i-1] + (dt/6.)*(k1[i-1] + 2*k2[i-1] + 2*k3[i-1] + k4[i-1])
        h[i] = h[i-1] + (dt/6.)*(l1[i-1] + 2*l2[i-1] + 2*l3[i-1] + l4[i-1])  
     
     # plot the solutions of the finite difference schemes 
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t*2, T*7.5, color='black', linewidth=2)  
    plt.plot(t*2, h*15, color = 'red', linewidth=1.5)
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,120])
    plt.ylim([-30,30])
    plt.legend()
    plt.xlabel('$Time(month)$')
    plt.ylabel('$T,h$')
    plt.savefig('plots/' + name_fig + '_time series.pdf')

    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(T*7.5, h*15, color='black', linewidth=2)  
    plt.axhline(0, linestyle=':', color='black')
    plt.axvline(0, linestyle=':', color='black')
    plt.xlim([-20,20])
    plt.legend()
    plt.xlabel('$T$')
    plt.ylabel('$h$')
    plt.savefig('plots/' + name_fig + '_trajactory.pdf')
    
    
def taskb(r = 0.25, alpha = 0.125, b = 2.5, N = 121, dt = 0.5, gama = 0.75, \
          miu, name_fig = 'TASK B'):

    # Set the variable lists
    T = np.zeros(N)
    h = np.zeros(N)
    t = np.zeros(N)
    # Define the Runge-Kutta scheme 
    k1 = np.zeros(N-1)
    k2 = np.zeros(N-1)
    k3 = np.zeros(N-1)
    k4 = np.zeros(N-1)
    l1 = np.zeros(N-1)
    l2 = np.zeros(N-1)
    l3 = np.zeros(N-1)
    l4 = np.zeros(N-1)

    T[0] = 1.125
    h[0] = 0.
    t[0] = 0.

    for i in range(1,N):
        t[i] = 0. + dt*i
    
        k1[i-1] = r*T[i-1] + gama*h[i-1]
        l1[i-1] = -r*h[i-1] - alpha*b*miu*T[i-1]
        k2[i-1] = r*(T[i-1] + (dt/2)*k1[i-1]) + gama*(h[i-1] + (dt/2)*l1[i-1])
        l2[i-1] = -r*(h[i-1] + (dt/2)*l1[i-1]) - \
                  alpha*b*miu*(T[i-1] + (dt/2)*k1[i-1])
        k3[i-1] = r*(T[i-1] + (dt/2)*k2[i-1]) + gama*(h[i-1] + (dt/2)*l2[i-1])
        l3[i-1] = -r*(h[i-1] + (dt/2)*l2[i-1]) - \
                  alpha*b*miu*(T[i-1] + (dt/2)*k2[i-1])
        k4[i-1] = r*(T[i-1] + dt*k3[i-1]) + gama*(h[i-1] + dt*l3[i-1])
        l4[i-1] = -r*(h[i-1] + dt*l3[i-1]) - alpha*b*miu*(T[i-1] + dt*k3[i-1])
    
        T[i] = T[i-1] + (dt/6.)*(k1[i-1] + 2*k2[i-1] + 2*k3[i-1] + k4[i-1])
        h[i] = h[i-1] + (dt/6.)*(l1[i-1] + 2*l2[i-1] + 2*l3[i-1] + l4[i-1])  
     
     # plot the solutions of the finite difference schemes 
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(t*2, T*7.5, color='black', linewidth=2)  
    plt.plot(t*2, h*15, color = 'red', linewidth=1.5)
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,120])
    plt.ylim([-30,30])
    plt.legend()
    plt.xlabel('$Time(month)$')
    plt.ylabel('$T,h$')
    plt.savefig('plots/' + name_fig + '_time series.pdf')

    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(T*7.5, h*15, color='black', linewidth=2)  
    plt.axhline(0, linestyle=':', color='black')
    plt.axvline(0, linestyle=':', color='black')
    plt.xlim([-20,20])
    plt.legend()
    plt.xlabel('$T$')
    plt.ylabel('$h$')
    plt.savefig('plots/' + name_fig + '_trajactory.pdf')