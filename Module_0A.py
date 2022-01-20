# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 15:00:28 2022

@author: Allison
"""

import numpy as np

# %%parameters with unit conversions
## INPUT the experimental system parameters
parameters = []
mu = float(input('Enter fluid viscosity (dyne-s/cm\u00b2): ')) * 1e-13
a = float(input('Enter cell/sphere radius (\u03BCm): ')) * 1e6
d = float(input('Enter critical distance (\u03BCm): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (\u03BCm): ')) * 1e6
b /= 2
w = float(input('Enter flow chamber width (\u03BCm): ')) * 1e6
parameters.append([mu, a, b, L, w, d])

## SET the experimental system parameters to reduce redundancy
#a = 7.85 * 1e6 # cell radius (um)
#d = 0.050 * 1e6 # distance from wall to edge of cell (um)
#L = 36 * 1e3 # bond length (nm)
#w = 800 * 1e6# c hamber width (um)
#b = 37.5 * 1e6# chamber half height (um)
#mu = 6.92e-3 * 1e-13# (dyne/cm^2)


Q = input('Enter flow rates (\u03BCL/hr), separated by commas and without new lines: ')
Q_str = [val.strip() for val in Q.split(',')]

Q_arr = np.zeros(len(Q_str))
for i in range(len(Q_str)):
    # Q_arr[i] = float(Q_str[i])  
    
    # convert from microliter/h to pm^3/s
    Q_arr[i] = float(Q_str[i]) * (10**27) / 3600 

Q_arr_nc = np.zeros(len(Q_str))
for i in range(len(Q_str)):
    Q_arr_nc[i] = float(Q_str[i])
    
# tether force
forces = Q_arr * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)

print('')
for i in range(len(Q_arr_nc)):
    print('For Q = %f (\u03BCL/hr), force = %f pN' % (Q_arr_nc[i], forces[i]))