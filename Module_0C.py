# -*- coding: utf-8 -*-
"""
Created on Mon May 31 11:54:40 2021

@author: Allison
"""

import pandas as pd
import numpy as np
import math
from scipy import stats
from matplotlib import pyplot as plt, cm
import os

#user directory/pathname
os.chdir('/User/RLNEK_tests/Module_0C_test')

# %% user input of data (single-molecule criteria)
#Compare m_l/C_l bond lifetimes at different flow rates (ONLY works for comapring 2 values of m_l/C_l)
files = [] #  spots files from Trackmate
m_l_list = [] # ASSUMES only 2 m_l values per Q (lists within lists)
C_l_list = [] # ASSUMES only 2 C_l values per Q (lists within lists)
Q_list = [] 


while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        # flow rate
        Q = float(input('Enter flow rate (\u03BCL/hr): '))
        Q_list.append(Q)

        # site densities or coating conc?
        which_one = input('Enter \"c\" to input coating concentrations or \"m\" to enter site densities: ')
        
        # characterized site densities
        if which_one == 'm':
            m_l_sublist = []
            m_l = input('For flow rate = %.2f (\u03BCL/hr), enter 2 site densities (sites/\u03BCm\u00b2) separated by commas and without new lines: ' % Q)
            m_l_str = [val.strip() for val in m_l.split(',')]
            
            for i in range(len(m_l_str)):
                m_l_sublist.append(float(m_l_str[i]))
                
            m_l_list.append(m_l_sublist)
            
            # spots in track statistics file per site density
            files_sublist = []
            for i in range(len(m_l_sublist)):
                spots_file = input('For flow rate = %.2f (\u03BCL/hr) and site density = %f (sites/\u03BCm\u00b2), enter name of "spots" file(s) from Trackmate: ' % (Q, m_l_sublist[i]))
                spots_file_list = [val.strip() for val in spots_file.split(',')]
                spots_file_subsub = []
                for j in range(len(spots_file_list)):
                    if not '.csv' in spots_file_list[j]:
                        spots_file_list[j] += '.csv'
                        
                    try:
                        with open(spots_file_list[j], encoding="unicode_escape") as file_open:
                            file = file_open.read()
                            spots_file_subsub.append(spots_file_list[j])
                    
                    except FileNotFoundError:
                        print('Invalid file name.')
                files_sublist.append(spots_file_subsub)
            files.append(files_sublist)
        
        elif which_one == 'c':
            # coating conc
            C_l_sublist = []
            C_l = input('For flow rate = %.2f (\u03BCL/hr), enter 2 coating concentrations (\u03BCg/mL) separated by commas and without new lines: ' % Q)
            C_l_str = [val.strip() for val in C_l.split(',')]
            
            for i in range(len(C_l_str)):
                C_l_sublist.append(float(C_l_str[i]))
            C_l_list.append(C_l_sublist)
        
            # 1 Spots in track statistics file for each C_l
            files_sublist = []
            for i in range(len(C_l_sublist)):
                spots_file = input('For flow rate = %.2f (\u03BCL/hr) and coating concentration = %f (\u03BCg/mL), enter name of "spots" file(s) from Trackmate: ' % (Q, C_l_sublist[i]))
                spots_file_list = [val.strip() for val in spots_file.split(',')]
                spots_file_subsub = []
                for j in range(len(spots_file_list)):
                    if not '.csv' in spots_file_list[j]:
                        spots_file_list[j] += '.csv'
                        
                    try:
                        with open(spots_file_list[j], encoding="unicode_escape") as file_open:
                            file = file_open.read()
                            spots_file_subsub.append(spots_file_list[j])
                    
                    except FileNotFoundError:
                        print('Invalid file name.')
                        
                    #spots_file_subsub.append(spots_file_list[i])
                
                # if not '.csv' in spots_file:
                #     spots_file += '.csv'
                
                files_sublist.append(spots_file_subsub)
                
            # files_sublist.append(spots_file_subsub)
            
            files.append(files_sublist)
                    
        else:
            print('Please enter \"c\" or \"m\".')
    
    else:
        print('Please type y/n.')
        
# %%user input of experimental system parameters

## INPUT the experimental system parameters
mu = float(input('Enter fluid viscosity (dyne-s/cm\u00b2): ')) * 1e-13
a = float(input('Enter cell/sphere radius (\u03BCm): ')) * 1e6
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

t_min_input = input('Enter non-specific binding time (enter for 0.2 seconds): t_min (seconds) = ')
if t_min_input == '':
    t_min_input = float(0.2)
else:
    t_min_input = float(t_min_input)
    
f_list = np.array(Q_list) * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)

# %% calculating bond lifetimes
def calc_disp(x0,x,y0,y):
        return np.sqrt((x-x0)**2+(y-y0)**2)
    
bond_times = []

for q in range(len(files)):
    
    lifetimes_per_Q = [] # lifetimes for each Q
    
    for n in range(len(files[q])):
        
        lifetimes_subsub = []
        
        for r in range(len(files[q][n])):
            spots_raw_data = pd.read_csv(files[q][n][r], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
    
            
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
            tmin_frames = math.ceil(t_min_input * CCD_FPS)
            while i < i_max-1:
                disp1 = calc_disp(x_pos[i+1], x_pos[j], y_pos[i+1], y_pos[j])
                if disp1 <= stop_dist:
                    i += 1
                    disp2 = calc_disp(x_pos[i], x_pos[j], y_pos[i], y_pos[j])
                    if i-j > tmin_frames:
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
            lifetimes_per_ml = []
            tc_particleID_new = []
            tc_trackID_new = []
            t_tot = 0
            
            # doing time conversion
            while i < len(tc_trackID):
                if tc_trackID[i] == tc_trackID[j]:
                    if tc_frame[i]-tc_frame[k] == 1:
                        t_tot += (tc_frame[i]-tc_frame[k]+6) / CCD_FPS
                        if i == len(tc_trackID)-1:
                            lifetimes_per_ml.append(t_tot)
                            tc_particleID_new.append(tc_particleID[k])
                            tc_trackID_new.append(tc_trackID[k])
                        i += 1
                        k += 1
                    else:
                        lifetimes_per_ml.append(t_tot)
                        tc_particleID_new.append(tc_particleID[k])
                        tc_trackID_new.append(tc_trackID[k])
                        t_tot = 0
                        j = i
                        i += 1
                        k += 1
                else:
                    lifetimes_per_ml.append(t_tot)
                    tc_particleID_new.append(tc_particleID[k])
                    tc_trackID_new.append(tc_trackID[k])
                    t_tot = 0
                    j = i
                    i += 1
                    k += 1  
            
            lifetimes_subsub.append(lifetimes_per_ml) 
            
        avg_lifetimes = []
        avg_lifetimes.append(np.mean(list(np.concatenate(lifetimes_subsub).flat)))
        lifetimes = []
        lifetimes.append(list(np.concatenate(lifetimes_subsub).flat))            
        lifetimes_per_Q.append(lifetimes)
        
    bond_times.append(lifetimes_per_Q)
    
# %% comparing lifetimes for 2 site densities/coating concs (for each Q)
lifetime_bin_vals = []  

threshold = 0.05 # for p < 0.05
t_stats = []
p_values = [] 

interval = input('Enter bin size for bond lifetimes (enter for 0.5 seconds): ')
if interval  == '':
    interval  = float(0.5)
else:
    interval  = float(interval)


# formatting plot
colors = iter(cm.rainbow(np.linspace(0, 1, len(Q_list))))
marks =['o', '+']
shapes = [marks for i in range(len(Q_list))] # different site densities

plt.figure(0)
plt.xlabel('Bond lifetimes (s)')
plt.ylabel('Fraction of stopping events')

for q in range(len(files)):

    lifetime_bin_vals_sublist = []
    c = next(colors)
    for f in range(len(files[q])):
        # bin lifetimes (first m_l or C_l)
        lifetime_series = pd.Series(list(np.concatenate(bond_times[q][f]).flat))
        num_bins = int(np.ceil((max(lifetime_series)) / interval))
        bins = lifetime_series.value_counts(normalize=True,sort=False,bins=num_bins)
        lifetime_bin_vals_sublist.append(bins.values)
            
        # plotting the distribution
        lifetimes_plt = np.linspace(0, max(lifetime_series), len(bins.values))
        
        
        # coating concentrations
        if len(C_l_list) != 0:
            plt.scatter(lifetimes_plt, bins.values, color=c, marker=shapes[q][f],
                     label='$Q = %.1f, C_l = %d$' % (Q_list[q], C_l_list[q][f]))
        else:
            plt.scatter(lifetimes_plt, bins.values, color=c, marker=shapes[q][f],
                     label='$Q = %.1f, m_l = %d$' % (Q_list[q], m_l_list[q][f]))

       
    # perform Welch's t-test
    # ASSUMES there will only be 2 site densities/coating conc per Q
    t_stat, pvalue = stats.ttest_ind(list(np.concatenate(bond_times[q][0]).flat),
                                     list(np.concatenate(bond_times[q][1]).flat),
                                     equal_var=False)
    t_stats.append(t_stat)
    p_values.append(pvalue)
    if len(m_l_list)>0:
        if (np.mean(list(np.concatenate(bond_times[q][0]).flat) >= np.mean(list(np.concatenate(bond_times[q][1]).flat)))) and (pvalue<threshold):
            print('')
            print('Bond lifetimes are significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in the site density range: ', end='')
            print(m_l_list[q], end='')
            print(u' (sites/\u03BCm\u00b2).')
            print('WARNING: This implies interactions may NOT meet the single-molecule criteria. Reassess input file(s)!')
            
        elif (np.mean(list(np.concatenate(bond_times[q][0]).flat) < np.mean(list(np.concatenate(bond_times[q][1]).flat)))) and (pvalue<threshold):
            print('')
            print('Bond lifetimes are significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in the site density range: ', end='')
            print(m_l_list[q], end='')
            print(u' (sites/\u03BCm\u00b2).')
            print('WARNING: This implies interactions may NOT meet the single-molecule criteria. Reassess input file(s)!')
            
        else:
            print('')
            print('Bond lifetimes are NOT significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in the site density range: ', end='')
            print(m_l_list[q], end='')
            print(u' (sites/\u03BCm\u00b2).') 
            print('This implies interactions may meet the single-molecule criteria in this range!')
            
    elif len(C_l_list)>0:
        if (np.mean(list(np.concatenate(bond_times[q][0]).flat) >= np.mean(list(np.concatenate(bond_times[q][1]).flat)))) and (pvalue<threshold):
            print('')
            print('Bond lifetimes are significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in coating concentration range: ', end='')
            print(C_l_list[q], end='')
            print(u' (\u03BCg/mL).')
            print('WARNING: This implies interactions may NOT meet the single-molecule criteria. Reassess input file(s)!')
            
        elif (np.mean(list(np.concatenate(bond_times[q][0]).flat) < np.mean(list(np.concatenate(bond_times[q][1]).flat)))) and (pvalue<threshold):
            print('')
            print('Bond lifetimes are significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in coating concentration range: ', end='')
            print(C_l_list[q], end='')
            print(u' (\u03BCg/mL).')
            print('WARNING: This implies interactions may NOT meet the single-molecule criteria. Reassess input file(s)!')
            
        else:
            print('')
            print('Bond lifetimes are NOT significantly different, ', end='')
            print('p=%f'% pvalue, end='')
            print(', for flow rate %f' %Q_list[q], end='')
            print(' (\u03BCL/hr) ', end='')
            print('in coating concentration range: ', end='')
            print(C_l_list[q], end='')
            print(u' (\u03BCg/mL).') 
            print('This implies interactions may meet the single-molecule criteria in this range!')
    else:
        print('Data was incorrectly input. Reassess input file(s)!')

plt.legend()
plt.savefig('single_molecule_criteria.png', dpi=300, bbox_inches='tight')