# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 09:33:49 2021

@author: Allison
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('/User/RLNEK_tests')

blue_circle_file = pd.read_csv('sm_blue_circles.csv')
blue_circle_time = blue_circle_file['Time (s)']
blue_circle_frac = blue_circle_file['Fraction']

blue_cross_file = pd.read_csv('sm_blue_crosses.csv')
blue_cross_time = blue_cross_file['Time (s)']
blue_cross_frac = blue_cross_file['Fraction']

red_circle_file = pd.read_csv('sm_red_circles.csv')
red_circle_time = red_circle_file['Time (s)']
red_circle_frac = red_circle_file['Fraction']

red_cross_file = pd.read_csv('sm_red_crosses.csv')
red_cross_time = red_cross_file['Time (s)']
red_cross_frac = red_cross_file['Fraction']

plt.figure(0)
plt.annotate(r'$Q_2$', 
             xy = (blue_circle_time[23], blue_circle_frac[23]),
             xytext=(0.65, 0.37), 
             color='b', 
             fontsize=19)
plt.plot(blue_circle_time, blue_circle_frac, 'bo', label=r'$m_{l,\,1}$')
plt.plot(blue_cross_time, blue_cross_frac, 'b+', label=r'$m_{l,\,2}$')

# plt.plot(purple_circle_time, purple_circle_frac, 'mo')
# plt.plot(purple_cross_time, purple_cross_frac, 'm+')

plt.annotate(r'$Q_1$', 
             xy = (red_circle_time[40], red_circle_frac[40]),
             xytext=(2.4, 0.72), 
             color='r', 
             fontsize=19)
plt.plot(red_circle_time, red_circle_frac, 'ro', label=r'$m_{l,\,1}$')
plt.plot(red_cross_time, red_cross_frac, 'r+', label=r'$m_{l,\,2}$')

plt.xlabel('Bond lifetimes (s)', fontsize=19)
plt.ylabel('Fraction of stopping events', fontsize=19)

ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])

plt.legend(fontsize=14)
plt.savefig('single_molecule.png', dpi=300, bbox_inches='tight')