# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 12:18:24 2021

@author: Allison
"""

import numpy as np
import scipy.optimize as optimize
import pandas as pd
import os

k_b = 0.0138
temp = 310.15 # Kelvin

os.chdir('/User/RLNEK_tests')

# %%
# functions
def catch_slip(f,E_21,k_1rup,k_2rup,f_12,x_B):
    exp_1 = np.exp(E_21/(k_b*temp))
    exp_2 = np.exp(f/f_12)
    exp_3 = np.exp((x_B*f)/(k_b*temp))
    return (exp_1*k_1rup + exp_2*k_2rup*exp_3) / (exp_1 + exp_2)

def slip(x,y):
    slope, log_koff_0 = np.polyfit(x,y,1)
    x_B = slope*k_b*temp
    koff_0 = np.exp(log_koff_0)
    return x_B, koff_0

# slip model k_off
def slip_func(x,x_B,koff_0):
    return koff_0*np.exp((x_B*x)/(k_b*temp))

# %%
# OVA ligand (red data from Liu et al, Figure 2B)
liu_red = pd.read_csv('liu_red.csv')
red_force = liu_red['Force']
red_time = 1 / liu_red['Time']

# catch-slip
red_params, red_matrix = optimize.curve_fit(catch_slip, red_force, red_time,
                                    bounds = np.array([0, np.inf]))

# slip 
red_x_B_s, red_koff_0 = slip(red_force, np.log(red_time))

# %%
# A2 ligand (blue data from Liu et al., Figure 2B)
liu_blue = pd.read_csv('liu_blue.csv')
blue_force = liu_blue['Force']
blue_time = 1 / liu_blue['Time']

# catch-slip
blue_params, blue_matrix = optimize.curve_fit(catch_slip, blue_force, blue_time,
                                    bounds = np.array([0, np.inf]))

# slip
blue_x_B_s, blue_koff_0 = slip(blue_force, np.log(blue_time))

# %%
# G4 ligand (purple data from Liu et al., Figure 2B)
liu_purple = pd.read_csv('liu_purple.csv')
purple_force = liu_purple['Force']
purple_time = 1 / liu_purple['Time']

# catch-slip
purple_params, purple_matrix = optimize.curve_fit(catch_slip, purple_force, purple_time,
                                                  bounds = np.array([0, np.inf]))

# slip
purple_x_B_s, purple_koff_0 = slip(purple_force, np.log(purple_time))

# %%
# R4 ligand (dark blue data from Liu et al., Figure 2C)
liu_db = pd.read_csv('liu_dark_blue.csv')
db_force = liu_db['Force']
db_time = 1 / liu_db['Time']

# catch-slip
db_params, matrix = optimize.curve_fit(catch_slip, db_force, db_time,
                                    bounds = np.array([0, np.inf]))

# slip
db_x_B_s, db_koff_0 = slip(db_force, np.log(db_time))

# %%
# E1 ligand (green data from Liu et al., Figure 2C)
liu_green = pd.read_csv('liu_green.csv')
green_force = liu_green['Force']
green_time = 1 / liu_green['Time']

# catch-slip
green_params, matrix = optimize.curve_fit(catch_slip, green_force, green_time,
                                    bounds = np.array([0, np.inf]))

# slip
green_x_B_s, green_koff_0 = slip(green_force, np.log(green_time))