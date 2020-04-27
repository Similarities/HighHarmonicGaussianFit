# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:11:46 2020

@author: similarities
"""
import numpy as np
from numpy import loadtxt
import math

from lmfit.models import GaussianModel

data = loadtxt('20190123_12 z1600 GVD900 PPinlineout.txt')
x = data[:, 0]
x[::] = x[::]-17.5/2
y = data[:, 1]
plt.figure(1)
plt.plot(x,y, label = 'data')



mod = GaussianModel()

pars = mod.guess(y, x=x)
out = mod.fit(y, pars, x=x)

print(out.fit_report(min_correl=0.15))

sigma = out.params['sigma'].value
print(sigma, type(sigma))

xx = np.linspace(-17.5/2,17.5/2, 1000)
yy = np.linspace(0,17.5, 1000)
#print(xx)

#/(2*(2.674**2)
for x in range(0,len(xx)): 
    yy[x] = (866363.685 /(2.674*((2*math.pi)**0.5)))* math.exp((-(xx[x]+2.25)**2)/(2*2.674**2))
    
plt.figure(1)
#plt.plot(x,y, label = 'data')
plt.plot(xx,yy, label = 'fit')
plt.legend()
plt.show()
    



