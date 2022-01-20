# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

@author: Allison
"""
import pandas as pd
import numpy as np
from scipy import stats
import os

os.chdir('/User/RLNEK_tests/Module_0B_test')
# %%

## INPUT the experimental system parameters
mu = float(input('Enter viscosity (dyne-s/cm\u00b2): ')) * 1e-13
a = float(input('Enter cell radius (\u03BCm): ')) * 1e6
d = float(input('Enter critical distance (\u03BCm): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (\u03BCm): ')) * 1e6
b /= 2
w = float(input('Enter flow chamber width (\u03BCm): ')) * 1e6
# parameters.append([mu, a, b, L, w, d])

## SET the experimental system parameters to reduce redundancy
#a = 7.85 * 1e6 # cell radius (um)
#d = 0.050 * 1e6 # distance from wall to edge of cell (um)
#L = 36 * 1e3 # bond length (nm)
#w = 800 * 1e6# c hamber width (um)
#b = 37.5 * 1e6# chamber half height (um)
#mu = 6.92e-3 * 1e-13# (dyne/cm^2)

y = a+d

CCD_FPS = float(input('Enter CCD FPS: '))
stop_dist = input(u'Enter minimum displacement (enter for 0.5 \u03BCm): D_min (\u03BCm) = ')
if stop_dist == '':
    stop_dist = float(0.5)
else:
    stop_dist = float(stop_dist)


# coating concs list
C_l_vals = []
#site density lists
spec_site_density = []
nonspec_site_density = []

# applied tensile/tether force list
spec_force = []
nonspec_force = []
idx = 1

# input coating concs/site densities 
check = input('Enter \"y\" if ligand site density was characterized; otherwise, enter \"n\": ')
if check.lower() == 'y':
    site_density = float(input('Enter ligand site density (sites/\u03BCm\u00b2): '))
    spec_site_density.append(site_density)
    nonspec_site_density.append(site_density)
    
elif check.lower() == 'n':
    C_l = float(input('Enter ligand coating concentration (\u03BCg/mL): '))
    C_l_vals.append(C_l)
    
# input flow rate
Q = float(input('Enter flow rate (\u03BCL/hr): '))

# input specific/nonspecific data    
nonspec_spots_data = input('For flow rate = %f (\u03BCL/hr), enter name of \"spots\" file(s) (from Trackmate) for non-specific ligand: ' % Q)
nonspec_spots_list = [val.strip() for val in nonspec_spots_data.split(',')]
nonspec_spots_files = []
for i in range(len(nonspec_spots_list)):
    if not '.csv' in nonspec_spots_list[i]:
        nonspec_spots_list[i] += '.csv'
        
    try:
        with open(nonspec_spots_list[i], encoding="unicode_escape") as file_open:
            file = file_open.read()
            nonspec_spots_files.append(nonspec_spots_list[i])
    
    except FileNotFoundError:
        print('Invalid file name.')

spec_spots_data = input('For flow rate = %f (\u03BCL/hr), enter name of \"spots\" file(s) (from Trackmate) for specific ligand: ' % Q)
spec_spots_list = [val.strip() for val in spec_spots_data.split(',')]
spec_spots_files = []
for i in range(len(spec_spots_list)):
    if not '.csv' in spec_spots_list[i]:
        spec_spots_list[i] += '.csv'
        
    try:
        with open(spec_spots_list[i], encoding="unicode_escape") as file_open:
            file = file_open.read()
            spec_spots_files.append(spec_spots_list[i])
    
    except FileNotFoundError:
        print('Invalid file name.')

# calculate applied tensile/tether force
Q = Q * (10**27) / 3600 
f = Q * np.sqrt(a/(2*L)) * (1.7005 * 9*np.pi*mu*a**2 + 0.9440 * 6*np.pi*mu*a**2) / (w*b**2)
spec_force.append(f)
nonspec_force.append(f)

# %% finding non-specific lifetime
# calculate displacement of a given cell/sphere
def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)

for i in range(len(nonspec_spots_files)):
    
    tmin_sublist = []
    spots_raw_data = pd.read_csv(nonspec_spots_files[i], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
            
    # r refers to meeting criteria
    r_pos_x = []
    r_pos_y = []
    r_trackID = []
    r_particleID = []
    r_frame = []
    
    i = 0
    j = 0
    
    x_pos = spots_raw_data['POSITION_X']
    y_pos = spots_raw_data['POSITION_Y']
    particle_ID = spots_raw_data['ID']
    track_ID = spots_raw_data['TRACK_ID']
    frames = spots_raw_data['FRAME']
    
    # number of iterations to calculate displacement
    i_max = len(frames) 
    
    # filter using stopping criteria
    while i < i_max-1:
        disp1 = calc_disp(x_pos[i+1], x_pos[j], y_pos[i+1], y_pos[j])
        if disp1 <= stop_dist:
            i += 1
            disp2 = calc_disp(x_pos[i], x_pos[j], y_pos[i], y_pos[j])
            if i-j > 6:
                r_particleID.append(particle_ID[i])
                r_trackID.append(track_ID[i])
                r_pos_x.append(x_pos[i])
                r_pos_y.append(y_pos[i])
                r_frame.append(frames[i])
        else:
            i += 1
            j = i-1
    
    # time conversion: (# of frames) -> seconds
    # tc = time conversion
    tc_particleID = np.array(r_particleID)
    tc_trackID = np.array(r_trackID)
    tc_frame = np.array(r_frame)
    
    # initial parameters
    i = 1
    j = 0
    k = 0
    nonspec_bond_times = []
    tc_particleID_new = []
    tc_trackID_new = []
    t_tot = 0
    
    # doing time conversion
    while i < len(tc_trackID):
        if tc_trackID[i] == tc_trackID[j]:
            if tc_frame[i]-tc_frame[k] == 1:
                t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                if i == len(tc_trackID)-1:
                    nonspec_bond_times.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                i += 1
                k += 1
            else:
                nonspec_bond_times.append(t_tot)
                tc_particleID_new.append(tc_particleID[k])
                tc_trackID_new.append(tc_trackID[k])
                t_tot = 0
                j = i
                i += 1
                k += 1
        else:
            nonspec_bond_times.append(t_tot)
            tc_particleID_new.append(tc_particleID[k])
            tc_trackID_new.append(tc_trackID[k])
            t_tot = 0
            j = i
            i += 1
            k += 1     
        
    t_min_avg = np.mean(nonspec_bond_times) # trial AVG
    tmin_sublist.append(nonspec_bond_times) # append list of bond lifetimes for given trial to master list

tot_nonspec_bond_times=list(np.concatenate(tmin_sublist).flat) # flatten list over trials
nonspec_tmin = np.mean(list(np.concatenate(tmin_sublist).flat)) # calculate mean over trials
nonspec_sem = stats.sem(list(np.concatenate(tmin_sublist).flat)) # calculate sem over trials
print('Non-specific bond lifetime (seconds): %f' % nonspec_tmin, end='')
print(u' \u00b1 %f' % nonspec_sem)

# finding specific lifetime
# track_raw_data = pd.read_csv(specific_track_data)
for i in range(len(spec_spots_files)):
    
    tmin_sublist = []
    spots_raw_data = pd.read_csv(spec_spots_files[i], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
    
    
    
    # r refers to meeting criteria
    # find which particles meet stopping criteria
    r_pos_x = []
    r_pos_y = []
    r_trackID = []
    r_particleID = []
    r_frame = []
    
    i = 0
    j = 0
    
    x_pos = spots_raw_data['POSITION_X']
    y_pos = spots_raw_data['POSITION_Y']
    particle_ID = spots_raw_data['ID']
    track_ID = spots_raw_data['TRACK_ID']
    frames = spots_raw_data['FRAME']
    
    # number of iterations to calculate displacement
    i_max = len(frames) 
    
    # filter using stopping criteria
    while i < i_max-1:
        disp1 = calc_disp(x_pos[i+1], x_pos[j], y_pos[i+1], y_pos[j])
        if disp1 <= stop_dist:
            i += 1
            disp2 = calc_disp(x_pos[i], x_pos[j], y_pos[i], y_pos[j])
            if i-j > 6:
                r_particleID.append(particle_ID[i])
                r_trackID.append(track_ID[i])
                r_pos_x.append(x_pos[i])
                r_pos_y.append(y_pos[i])
                r_frame.append(frames[i])
        else:
            i += 1
            j = i-1
    
    # time conversion: (# of frames) -> seconds
    # tc = time conversion
    tc_particleID = np.array(r_particleID)
    tc_trackID = np.array(r_trackID)
    tc_frame = np.array(r_frame)
    
    # initial parameters
    i = 1
    j = 0
    k = 0
    spec_bond_times = []
    tc_particleID_new = []
    tc_trackID_new = []
    t_tot = 0
    
    # doing time conversion
    while i < len(tc_trackID):
        if tc_trackID[i] == tc_trackID[j]:
            if tc_frame[i]-tc_frame[k] == 1:
                t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                if i == len(tc_trackID)-1:
                    spec_bond_times.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                i += 1
                k += 1
            else:
                spec_bond_times.append(t_tot)
                tc_particleID_new.append(tc_particleID[k])
                tc_trackID_new.append(tc_trackID[k])
                t_tot = 0
                j = i
                i += 1
                k += 1
        else:
            spec_bond_times.append(t_tot)
            tc_particleID_new.append(tc_particleID[k])
            tc_trackID_new.append(tc_trackID[k])
            t_tot = 0
            j = i
            i += 1
            k += 1
            
    t_min_avg = np.mean(spec_bond_times)#trial AVG
    tmin_sublist.append(spec_bond_times)#append list of bond lifetimes for given trial to master list
    
tot_spec_bond_times = list(np.concatenate(tmin_sublist).flat)#flatten list over trials
spec_tmin = np.mean(list(np.concatenate(tmin_sublist).flat))#calculate mean over trials
spec_sem = stats.sem(list(np.concatenate(tmin_sublist).flat))#calculate sem over trials
print('Specific bond lifetime (seconds): %f' % spec_tmin, end='')
print(u' \u00b1 %f' % spec_sem)

# %% Welch's t-test on specific and non-specific bond lifetimes
# should only be comparing 2 ligands @ equal flow rate and coating concentration/site density
t_stat, pvalue = stats.ttest_ind(tot_nonspec_bond_times,
                                  tot_spec_bond_times,
                                  equal_var=False)
if (nonspec_tmin >= spec_tmin) and (pvalue<0.05):
    print('p=%f'% pvalue)
    print('WARNING: Specific bond lifetimes are significantly LESS than non-specific bond lifetimes. Reassess input file(s)!')
    
elif (nonspec_tmin < spec_tmin) and (pvalue<0.05):
    print('p=%f'% pvalue)
    print('Specific bond lifetimes are significantly GREATER than non-specific bond lifetimes. Your data looks reasonable!')
    
else:
    print('p=%f'% pvalue)
    print('WARNING: Specific bond lifetimes are NOT significantly different than non-specific bond lifetimes. Reassess input file(s)!')

    
