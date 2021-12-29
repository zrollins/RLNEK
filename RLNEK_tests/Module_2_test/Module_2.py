#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 08:34:05 2020

@author: Allison
"""
import pandas as pd
import numpy as np
from scipy import stats, optimize, special, interpolate
from matplotlib import pyplot as plt, cm
import math
import os
import warnings
import csv
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 


os.chdir('/User/RLNEK_tests/Module_2_test')

# %%input experimental system parameters

parameters = []
##INPUT the experimental system parameters
mu = float(input('Enter fluid viscosity (dyne-s/cm\u00b2): ')) * 1e-13
a = float(input('Enter cell/sphere radius (\u03BCm): ')) * 1e6
d = float(input('Enter critical distance (\u03BCm): ')) * 1e6
L = float(input('Enter receptor-ligand bond length (nm): ')) * 1e3
b = float(input('Enter flow chamber height (\u03BCm): ')) * 1e6
b /= 2
w = float(input('Enter flow chamber width (\u03BCm): ')) * 1e6
parameters.append([mu, a, b, L, w, d])

##SET the experimental system parameters to reduce redundancy
#a = 7.85 * 1e6 # cell radius (um)
#d = 0.050 * 1e6 # distance from chamber floor to edge of cell (um)
#gamma_w = 1 # wall shear rate 
#L = 36 * 1e3 # bond length (nm)
#w = 800 * 1e6# c hamber width (um)
#b = 37.5 * 1e6# chamber half height (um)
#mu = 6.92e-3 * 1e-13# (dyne/cm^2)

y = a+d #critical disance of cell/sphere from chamber floor

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
    
m_r = float(input('Enter cell/sphere site density (sites/\u03BCm\u00b2): '))


# inter/extrapolation of the speed constant: hydrodynamic theoretic prediction of cell/sphere rotational to hydrodynamic velocity
da = np.array([0, 1e-8, 1e-7, 1e-6, 
                1e-5, 1e-4, 1e-3, 0.003202])
speed_constants = np.array([0.5676, 0.5556, 0.5539, 0.5518, 
                            0.5488, 0.5445, 0.5375, 0.5315])
fit_func = interpolate.interp1d(da, speed_constants,
                                fill_value='extrapolate')

speed_const = fit_func(d/a)

# %% loop that inputs user data
track_data = []
spots_data = []
site_densities = []
shear_rates = []
forces = []
Q_vals = []
N_T_vals = []
t_min_vals = []

while True:
    run = str(input('Enter \"y\" to input data or \"n\" to execute: ')) 
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        
        # site density
        m_l = float(input('Enter flow chamber site density (sites/\u03BCm\u00b2): '))
        site_densities.append(m_l)
        
        # flow rate
        Q = input('For flow chamber site density = %f (sites/\u03BCm\u00b2), enter flow rates (\u03BCL/hr), separated by commas and without new lines: ' % m_l)
        Q_str = [val.strip() for val in Q.split(',')]
        
        Q_arr = np.zeros(len(Q_str))
        for i in range(len(Q_str)):
            # Q_arr[i] = float(Q_str[i])  
            # convert from microliter/h to pm^3/s
            Q_arr[i] = float(Q_str[i]) * (10**27) / 3600 
        
        Q_arr_nc = np.zeros(len(Q_str))
        for i in range(len(Q_str)):
            Q_arr_nc[i] = float(Q_str[i])
        Q_vals.append(list(Q_arr_nc))
            
        # Calculate applied tensile/tether force
        f = Q_arr * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)
        forces.append(list(f))
        
        #Calculate shear stress & shear rate
        tau = (3*mu*Q_arr) / (2*w*b**2)
        shear_rate = tau / mu
        shear_rates.append(shear_rate)
        
        #Input Trackmate "tracks" and "spots" file(s) for conditions [site density][flow rate][trial]
        track_data_sublist = []
        spots_data_sublist = []
        t_min_sublist = []
        N_T_sublist = []
        
        for i in range(len(Q_arr)):
            track_file_name = input('For flow rate = %.2f (\u03BCL/hr) and flow chamber site density = %f (sites/\u03BCm\u00b2), enter name of "tracks" file(s) from Trackmate: ' % (Q_arr_nc[i], m_l))
            track_file_list = [val.strip() for val in track_file_name.split(',')]
            track_file_subsub = []
            for j in range(len(track_file_list)):
                if '.csv' not in track_file_list[j]:
                    track_file_list[j] += '.csv'
                
                try:
                    with open(track_file_list[j], encoding="unicode_escape") as file_open:
                        file = file_open.read()
                        track_file_subsub.append(track_file_list[j])
                
                except FileNotFoundError:
                    print('Invalid file name.')
                    
            track_data_sublist.append(track_file_subsub)
                
            spots_file_name = input('For flow rate = %.2f (\u03BCL/hr) and flow chamber site density = %f (sites/\u03BCm\u00b2), enter name of "spots" file(s) from Trackmate: ' % (Q_arr_nc[i], m_l))
            spots_file_list = [val.strip() for val in spots_file_name.split(',')]
            spots_file_subsub = []
            for j in range(len(spots_file_list)):
                if '.csv' not in spots_file_list[j]:
                    spots_file_list[j] += '.csv'
                    
                try:
                    with open(spots_file_list[j], encoding="unicode_escape") as file_open:
                        file = file_open.read()
                        spots_file_subsub.append(spots_file_list[j])
                
                except FileNotFoundError:
                    print('Invalid file name.')
                    
            spots_data_sublist.append(spots_file_subsub)
            
            N_T_input = input('For each trial of flow rate = %.2f (\u03BCL/hr) and flow chamber site density = %f (sites/\u03BCm\u00b2), enter N_T values(s): ' % (Q_arr_nc[i], m_l))
            N_T_str = [val.strip() for val in N_T_input.split(',')]
            N_T_subsub = []
            
            for j in range(len(N_T_str)):
                N_T_float = float(N_T_str[j])
                N_T_subsub.append(N_T_float)
                
            N_T_sublist.append(N_T_subsub)
        
        track_data.append(track_data_sublist)
        spots_data.append(spots_data_sublist)
        #t_min_vals.append(t_min_sublist)
        N_T_vals.append(N_T_sublist)

    else:
        print('Please enter \"y\" or \"n\".')

# %% calculating k_off
u_f_vals = []
U_cell_vals = []
U_hd_vals = []
U_cell_avg_vals = []
U_hd_avg_vals = []

koff_vals = []
koff_trackID_vals = []
koff_avg_vals = []
koff_error_vals = []

Nb_vals = []
NbNT_vals = []
NbNT_error_vals = []


#3-dimensional nested list format:
##sub = site density
##subsub = site density, flow rate
##subsubsub = site density, flow rate, trial
## e.g., Nb_vals[site density][flow rate][trial]
##AVG and SEM is over all inputted trials for a given condition [site density][flow rate]

for m in range(len(track_data)): 
    u_f_sub = []
    U_cell_sub = []
    U_hd_sub = []
    U_cell_avg_sub = []
    U_hd_avg_sub = []
    
    Nb_vals_sub = []
    NbNT_vals_sub = []
    NbNT_error_vals_sub = []
    
    koff_all_sub = []
    koff_trackID_vals_sub = []
    koff_avg_vals_sub = []
    koff_error_vals_sub = []
    
    for n in range(len(track_data[m])):
    #loop for multiple ligand site densities
        #critical velocity calculation for filtering tracks
        u_f = y*shear_rates[m][n]*(1-(5/16)*(a/y)**3) * 1e-6 # convert back to microns
        u_f_sub.append(u_f)
        U_cell_subsub = []
        U_hd_subsub = []

        koff_all_subsub = []
        koff_trackID_vals_subsub = []
        koff_avg_vals_subsub = []
        koff_error_vals_subsub = []
        
        Nb_vals_subsub = []
        NbNT_vals_subsub = []

        for p in range(len(track_data[m][n])):
        #loop for multiple trials
            #extract data into pandas
            tracks_raw_data = pd.read_csv(track_data[m][n][p], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
            spots_raw_data = pd.read_csv(spots_data[m][n][p], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
            
            #collect track velocities < critical velocity
            filtered_speeds = tracks_raw_data[tracks_raw_data['TRACK_MEAN_SPEED'] < np.absolute(u_f)]
            filtered_tracks_list = list(filtered_speeds['TRACK_ID'])
            
            #collect all trackIDs (with tracks < critical velocity) in spots file
            better_tracks = []
            trackID = spots_raw_data['TRACK_ID']
            particleID = spots_raw_data['ID']
            x_pos = spots_raw_data['POSITION_X']
            y_pos = spots_raw_data['POSITION_Y']
            frame = spots_raw_data['FRAME']
            
            for i in range(len(filtered_tracks_list)): 
                for j in range(len(trackID)):
                    if trackID[j] == filtered_tracks_list[i]:
                        if j != 0:
                            if trackID[j-1] != trackID[j]:
                                better_tracks.append(trackID[j])
                        else:
                            better_tracks.append(trackID[j])
                            
            #collect all particleIDs' data (with trackIDs < critical velocity) in spots file
            particleID_new = []
            trackID_new = []
            x_new = []
            y_new = []
            frame_new = []
            
            for i in range(len(better_tracks)):
                for j in range(len(trackID)):
                    if trackID[j] == better_tracks[i]:
                        particleID_new.append(particleID[j])
                        trackID_new.append(trackID[j])
                        x_new.append(x_pos[j])
                        y_new.append(y_pos[j])
                        frame_new.append(frame[j])
            
            #determine filtered particles that meet stopping criteria (D_min, t_min)
            #r refers to meeting criteria
           
            r_pos_x = []
            r_pos_y = []
            r_trackID = []
            r_particleID = []
            r_frame = []
            
            i = 1
            i_max = len(trackID_new)
            j = 0
            
            #calculate particle displacement, D
            def calc_disp(x0,x,y0,y):
                return np.sqrt((x-x0)**2+(y-y0)**2)
            
            #iterate through frames, calculate D, collect particle data that meet stopping criteria (D_min, t_min)
            tmin_frames = math.ceil(t_min_input * CCD_FPS)
            while i < i_max-1:
                disp1 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
                if disp1 <= stop_dist:
                    if i-j > tmin_frames:
                        r_particleID.append(particleID_new[i])
                        r_trackID.append(trackID_new[i])
                        r_pos_x.append(x_new[i])
                        r_pos_y.append(y_new[i])
                        r_frame.append(frame_new[i])
                    i += 1
                else:
                    i += 1
                    j = i-1
            
            #stopping events time conversion: (# of frames) -> seconds
            #tc = time conversion
            tc_particleID = np.array(r_particleID)
            tc_trackID = np.array(r_trackID)
            tc_frame = np.array(r_frame)
            tc_pos_x =np.array(r_pos_x)
            tc_pos_y =np.array(r_pos_y)
            
            #initial parameters
            t_total = []
            tc_trackID_new = []
            i = 1
            j = 0
            t_tot = 0
            
            #time conversion
            while i < len(tc_trackID):
                if tc_trackID[i] == tc_trackID[j]:
                    disp2 = calc_disp(tc_pos_x[i],tc_pos_x[j],tc_pos_y[i],tc_pos_y[j])
                    if ((tc_frame[i]-tc_frame[j] > 0) and (disp2 <= stop_dist)):
                        if  i == len(tc_trackID)-1:
                            t_tot += (tc_frame[i] - tc_frame[j])
                            t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                            tc_trackID_new.append(tc_trackID[j])
                            t_tot = 0
                            j=i  
                            i+=1
                        else:
                            t_tot += (tc_frame[i] - tc_frame[j])
                            j=i
                            i+=1
                    elif ((tc_frame[i]-tc_frame[j] > 0) and (disp2 > stop_dist)):
                        t_total.append((t_tot + 1) / CCD_FPS)
                        tc_trackID_new.append(tc_trackID[j])
                        t_tot = 0
                        j=i  
                        i+=1
                    else:
                        t_tot = 0
                        j=i
                        i+=1
                else:
                    t_total.append((t_tot + tmin_frames + 1) / CCD_FPS)
                    tc_trackID_new.append(tc_trackID[j])
                    t_tot=0
                    j=i
                    i +=1
                    
            #determine stopping events with unique track IDs, Nb
            i = 1
            j = 0
            k = 0
            t_total_unique = []
            tc_trackID_unique = []
            t_tot = np.array([0])
            
            while i < len(tc_trackID_new):
                if tc_trackID_new[i] != tc_trackID_new[j]:
                    tc_trackID_unique.append(tc_trackID_new[j])
                    t_tot = np.add(t_tot, t_total[k])
      
                    t_total_unique.append(t_tot)
                    if i == len(tc_trackID_new) - 1:
                        tc_trackID_unique.append(tc_trackID_new[i])
                        t_total_unique.append(t_total[i])
                        
                    j = i
                    i += 1
                    k += 1
                    t_tot = np.array([0])   
                else:
                    t_tot = np.add(t_tot, t_total[k])
                    
                    if i == len(tc_trackID_new) - 1:
                        t_tot = np.add(t_tot, t_total[i])
                        tc_trackID_unique.append(tc_trackID_new[i])
                        t_total_unique.append(t_tot)
   
                    i += 1
                    k += 1
            
            #calculate moving time of particles not bound, U_hd (for k+*)
            durations_raw = tracks_raw_data['TRACK_DURATION'] # units: seconds
            disp_raw = tracks_raw_data['TRACK_DISPLACEMENT'] # units: microns
            trackID_raw = tracks_raw_data['TRACK_ID']
            track_speeds_raw = tracks_raw_data['TRACK_MEAN_SPEED'] # units: microns/second
            
            #getting unique track information
            tracks_unique_data = tracks_raw_data.loc[(tracks_raw_data['TRACK_ID'].isin(tc_trackID_unique))]
            durations_unique = tracks_unique_data['TRACK_DURATION'].reset_index(drop=True)
            disp_unique2 = tracks_unique_data['TRACK_DISPLACEMENT'].reset_index(drop=True)
            track_speeds_unique =tracks_unique_data['TRACK_MEAN_SPEED'].reset_index(drop=True)
            
            time_moving_vals = []
            durations_sublist = []
            disp_unique = []
            U_hd_subsubsub = []
            U_cell_subsubsub = []
            for q in range(len(tc_trackID_unique)):
                U_cell = track_speeds_unique[q]
                durations_sublist.append(durations_unique[q])
                time_moving = durations_unique[q] - t_total_unique[q]
                time_moving_vals.append(time_moving)
                disp_unique.append(disp_unique2[q])
                U_hd = disp_unique2[q] / time_moving_vals[q]
                U_hd_subsubsub.append(U_hd)
                U_cell_subsubsub.append(U_cell)

            #update list of Uhd, Ucell per condition[site density][flow rate][trial]
            U_hd_subsubsub = [float(i) for i in U_hd_subsubsub]  #convert np.array to float  
            U_hd_subsub.append(U_hd_subsubsub)
            U_cell_subsub.append(U_cell_subsubsub)

            
            #compute k_off = 1/lifetimes
            k_off_subsubsub = []
            #k_off_subsubsub = np.reciprocal(t_total_unique, where=t_total_unique!=0)
            k_off_subsubsub = np.reciprocal(t_total, where=t_total!=0)
            k_off_subsubsub = [float(i) for i in k_off_subsubsub]
            
            #update list of all k_offs per condition [site density][flow rate][trial]
            koff_all_subsub.append(k_off_subsubsub)
            koff_trackID_vals_subsub.append(tc_trackID_new)

            #update list of tracks bound, N_b, per condition [site density][flow rate][trial]
            Nb = len(tc_trackID_unique) #number of tracks that meet stop criteria
            Nb_vals_subsub.append(Nb)
            
            #compute and update list of capture efficiency, Nb/NT, per condition [site density][flow rate][trial]
            NT = N_T_vals[m][n][p] #input NT values
            NbNT_vals_subsub.append(Nb/NT)  #nested list of Nb/NT values for single condition and n TRIALS
        
        #compute koff, Nb/NT AVG and SEM across all trials [site density][flow rate]
        koff_avg_new = np.mean(list(np.concatenate(koff_all_subsub).flat)) #flatten koff list across TRIALS and AVG for [site density][flow rate]
        koff_avg_error = stats.sem(list(np.concatenate(koff_all_subsub).flat)) # flatten koff list across TRIALS and SEM for [site density][flow rate]
        NbNT_avg_new = np.mean(NbNT_vals_subsub) # AVG Nb/NT across TRIALS for [site density][flow rate]
        NbNT_avg_error = stats.sem(NbNT_vals_subsub) #SEM Nb/NT across TRIALS for [site density][flow rate]
        
        #update lists of AVG and SEM 
        koff_avg_vals_sub.append(koff_avg_new) #append koff AVG for [site density][flow rate]
        koff_error_vals_sub.append(koff_avg_error) #append koff SEM for [site density][flow rate]
        NbNT_error_vals_sub.append(NbNT_avg_error) #append Nb/NT AVG for [site density][flow rate]
        NbNT_vals_sub.append(NbNT_avg_new) #append Nb/NT SEM for [site density][flow rate]
        U_hd_avg_sub.append(np.mean(list(np.concatenate(U_hd_subsub).flat)))  #append Uhd trial flattened AVG for [site density][flow rate]
        U_cell_avg_sub.append(np.mean(list(np.concatenate(U_cell_subsub).flat)))  #append Ucell trial flattened AVG for [site density][flow rate]
        
        #update all values lists
        koff_all_sub.append(koff_all_subsub)
        koff_trackID_vals_sub.append(koff_trackID_vals_subsub)
        Nb_vals_sub.append(Nb_vals_subsub)
        U_hd_sub.append(U_hd_subsub)
        U_cell_sub.append(U_cell_subsub)

    u_f_vals.append(u_f_sub)
    U_hd_vals.append(U_hd_sub)
    U_cell_vals.append(U_cell_sub)
    U_hd_avg_vals.append(U_hd_avg_sub)
    U_cell_avg_vals.append(U_cell_avg_sub)
    
    Nb_vals.append(Nb_vals_sub)
    NbNT_vals.append(NbNT_vals_sub)
    NbNT_error_vals.append(NbNT_error_vals_sub)
    
    koff_vals.append(koff_all_sub)
    koff_trackID_vals.append(koff_trackID_vals_sub)
    koff_avg_vals.append(koff_avg_vals_sub)
    koff_error_vals.append(koff_error_vals_sub)
    
# %% k_off fitting
k_b = 0.0138 #boltzmann constant in (J/K)e-21 or (kg*pm*nm)/(K*s^2)
temp = 310.15 # biological temperature in  Kelvin

#slip model k_off
def slip(x,y):
    slope, log_k_off_0 = np.polyfit(x,y,1)
    x_B = slope*k_b*temp
    k_off_0 = np.exp(log_k_off_0)
    return x_B, k_off_0

def slip_func(f,x_B,k_off_0):
    return k_off_0*np.exp((x_B*f)/(k_b*temp))

#catch-slip model k_off
def catch_slip(f,E_21,k_1rup,f_12,k_2rup,x_B):
    exp_1 = np.exp(E_21/(k_b*temp))
    exp_2 = np.exp(f/f_12)
    exp_3 = np.exp((x_B*f)/(k_b*temp))
    return (exp_1*k_1rup + exp_2*k_2rup*exp_3) / (exp_1 + exp_2)

#catch-slip fitting
#fit parameters: [E_21,k_1rup,f_12,k_2rup,x_B]
E_21_list = []
k_1rup_list = []
f_12_list = []
k_2rup_list = []
x_B_list_cs = []

for i in range(len(forces)):
    params, matrix = optimize.curve_fit(catch_slip, forces[i], koff_avg_vals[i], bounds=[0,np.inf])
    E_21_list.append(params[0])
    k_1rup_list.append(params[1])
    f_12_list.append(params[2])
    k_2rup_list.append(params[3])
    x_B_list_cs.append(params[4])

#slip fitting
x_B_list_s = []
koff_0_list = []

for i in range(len(forces)):
    x_B_s, koff_0 = slip(forces[i], np.log(koff_avg_vals[i]))
    #nonspec_slip_fit = slip_func(forces, x_B_nonspec, koff_0_nonspec)
    x_B_list_s.append(x_B_s)
    koff_0_list.append(koff_0)
    
#distinguishing correct bond dissociation model with physical constraint on x_B
def compare_x_B(x_B_s_arr, x_B_cs_arr):
    
    # catch-slip parameters
    x_B_cs_final = []
    E_21_final = []
    k_1rup_final = []
    k_2rup_final = []
    f_12_final = []
    
    # slip model parameters
    k_off_0_final = []
    x_B_slip_final = []
    
    s_forces = [] # force values (slip)
    cs_forces = [] # force values (catch-slip)
    s_k_off = [] # k_off values (slip)
    cs_k_off = [] # k_off values (catch-slip)

    if len(x_B_s_arr) != len(x_B_cs_arr):
        print('Arrays must be equal in length.')
        
    else:
        for i in range(len(x_B_cs_arr)):
            
            # if x_B is outside of a given threshold, the data is slip
            if (x_B_cs_arr[i] < 10**(-3)) or (x_B_cs_arr[i] > 10):
                x_B_slip_final.append(x_B_s_arr[i])
                k_off_0_final.append(koff_0_list[i])
                s_forces.append(forces[i])
                s_k_off.append(koff_avg_vals[i])
              
            # if x_B is inside a given threshold, the data is catch-slip
            else:
                x_B_cs_final.append(x_B_cs_arr[i])
                E_21_final.append(E_21_list[i])
                k_1rup_final.append(k_1rup_list[i])
                k_2rup_final.append(k_2rup_list[i])
                f_12_final.append(f_12_list[i])
                cs_forces.append(forces[i])
                cs_k_off.append(koff_avg_vals[i])
        
        # all slip
        if (len(x_B_slip_final) != 0) and (len(x_B_cs_final) == 0):
            return k_off_0_final, x_B_slip_final, s_forces, s_k_off
    
        # all catch-slip
        elif (len(x_B_slip_final) == 0) and (len(x_B_cs_final) != 0):
            return x_B_cs_final, E_21_final, k_1rup_final, k_2rup_final, f_12_final, cs_forces, cs_k_off
     
        elif (len(x_B_slip_final) != 0) and (len(x_B_cs_final) != 0):
            return x_B_cs_final, E_21_final, k_1rup_final, k_2rup_final, f_12_final, cs_forces, cs_k_off, k_off_0_final, x_B_slip_final, s_forces, s_k_off

#results of distinguishing
lists = compare_x_B(x_B_list_s, x_B_list_cs)
if len(lists) == 4:
    koff_0_vals = lists[0]
    x_B_vals = lists[1]
    force_vals_s = lists[2]
    koff_vals_s = lists[3]
    
elif len(lists) == 7:
    x_B_vals_cs = lists[0]
    E_21_vals = lists[1]
    k_1rup_vals = lists[2]
    k_2rup_vals = lists[3]
    f_12_vals = lists[4]
    force_vals_cs = lists[5]
    koff_vals_cs = lists[6]
    
elif len(lists) == 11:
    x_B_vals_cs = lists[0]
    E_21_vals = lists[1]
    k_1rup_vals = lists[2]
    k_2rup_vals = lists[3]
    f_12_vals = lists[4]
    force_vals_cs = lists[5]
    koff_vals_cs = lists[6]
    koff_0_vals = lists[7]
    x_B_vals_s = lists[8]
    force_vals_s = lists[9]
    koff_vals_s = lists[10]
    
#k+ function
def k_plus_func(NbNT, k_off, m_l):
    num = NbNT*k_off
    dem = m_l*(1-NbNT)
    return num/dem

#compute k+ over conditions [site density][flow rate]
kplus_vals = []
for i in range(len(NbNT_vals)):
    kplus = k_plus_func(np.array(NbNT_vals[i]),
                        koff_avg_vals[i],
                        site_densities[i])
    #error = stats.sem(kplus)
    kplus_vals.append(kplus)
   
#compute AVG, SEM k+ over conditions [site density]
kplus_vals_avg = []
kplus_errors = []
for i in range(len(kplus_vals[0])):
    avg_sublist = [sublist[i] for sublist in kplus_vals]
    kplus_vals_avg.append(np.mean(avg_sublist))
    kplus_errors.append(stats.sem(avg_sublist))
    
#compute k_in

#set system parameters
#D = 0.15 # receptor diffusivity (um^2/s)
#alpha = 2 # reactive radius (nm)

#input system parameters
D = float(input('Enter receptor-ligand diffusivity (\u03BCm\u00b2/s): '))  # um^2/s
alpha = float(input('Enter ligand reactive radius (nm): ')) * 1e-3 # convert to micron

#Peclet number function
def Pe_func(u_f_arr,alpha,D):
    u_f_arr *= np.asarray(1-speed_const) 
    Pe = u_f_arr*alpha/D
    return Pe

#Nusselt number function
def Nu_func(Pe):
    I_0 = special.iv(0,Pe/2)
    K_0 = special.kv(0,Pe/2)
    summation_Nu = 0

    for n in range(1,100):
        I_n = special.iv(n,Pe/2)
        K_n = special.kv(n,Pe/2)
        summation_Nu += (-1)**n*(I_n/K_n)
    Nu = 2*((I_0/K_0) + 2*summation_Nu)
    
    return Nu

#Hammer number function (dimensionless duration time)
def lambda_func(Pe):
    I_0 = special.iv(0,Pe/2)
    I_1 = special.iv(1,Pe/2)
    fraction_1 = -I_1**3/I_0
    summation_lambda = 0
    
    for n in range(1,74): # max: 74, otherwise: runtime warning
        I_smol = special.iv(n-1,Pe/2)
        I_large = special.iv(n+1,Pe/2)
        I_n = special.iv(n,Pe/2)
        num = I_smol*I_large*(I_smol+I_large)
        dem = I_n
        summation_lambda += ((-1)**(n+1))*(num/dem)

    cap_lambda = (1/Pe)*(fraction_1 + summation_lambda)
    
    return cap_lambda

#Damkohler number function
def k_plus_Pe(Pe,delta, Ac):
    Nu = Nu_func(Pe)
    cap_lambda = lambda_func(Pe)
    P = (cap_lambda*delta)/(1+cap_lambda*delta)
    return np.pi*m_r*Ac*D*Nu*P

#k_in function
def k_in_func(alpha,D,delta):
    return (delta*D)/alpha**2

#compute Pe numbers
Pe_vals = []
for i in range(len(u_f_vals)):
    Pe = []
    Pe = Pe_func(u_f_vals[i], alpha, D)
    Pe_vals.append(list(Pe))

#compute Nu numbers
Nu_vals = []
for i in range(len(Pe_vals)):
    Nu = Nu_func(np.array(Pe_vals[i]))
    Nu_vals.append(list(Nu))

#compute Ha (lambda) numbers    
lambda_vals = []    
for i in range(len(Pe_vals)):
    lamb = lambda_func(np.array(Pe_vals[i]))
    lambda_vals.append(list(lamb))


#calculate cell/sphere max contact area
Ac_max = 4*np.pi*(a*1e-6)**2
delta_vals = []
Ac_vals = []

#compute Da numbers
for i in range(len(Pe_vals)):
    params, matrix = optimize.curve_fit(k_plus_Pe, Pe_vals[i], kplus_vals[i], bounds = [[0,0],[np.inf,Ac_max]])
    delta_vals.append(params[0])
    Ac_vals.append(params[1])

#compute k_in for each condition [site density]   
kin_vals = []
for i in range(len(delta_vals)):
    kin = k_in_func(alpha, D, delta_vals[i])
    kin_vals.append(kin)

#AVG, SEM Ac over conditions [site density]
Ac_avg = np.mean(Ac_vals)
Ac_err = stats.sem(Ac_vals)

#AVG, SEM k_in over conditions [site density]     
kin_avg = np.mean(kin_vals)
kin_err = stats.sem(kin_vals)

# %% k+* validation

#k+* function
def kplus_star_func(U_hd,U_cell,koff,m_l):
    ratio = np.array(U_cell) / np.array(U_hd)
    kplus = (np.array(koff)/m_l) * ((1-ratio)/ratio)
    return kplus

#compute AVG, SEM k+* over conditions [site density][flow rate]
kplus_star_vals = []
for i in range(len(koff_avg_vals)):
    kplus_star = kplus_star_func(U_hd_avg_vals[i], 
                                 U_cell_avg_vals[i], 
                                 koff_avg_vals[i],
                                 site_densities[i])
    kplus_star_vals.append(kplus_star)


#compute AVG, SEM k+* over conditions [site density]
kplus_star_avg = []
kplus_star_errors = []
for i in range(len(kplus_star_vals[0])):
    avg_sublist = [sublist[i] for sublist in kplus_star_vals]
    kplus_star_avg.append(np.mean(avg_sublist))
    kplus_star_errors.append(stats.sem(avg_sublist))   
    

# %% curve fitting: k_off, k+, Nb/NT
f_fit_vals = np.linspace(forces[0][0], forces[0][-1], 1000)

#koff fit function
def koff_fit_func(f):
    koff_fit = []
    for i in range(len(site_densities)):
        # slip model
        if (len(lists) == 4) or (len(lists) == 11):
            koff_s = slip_func(f, x_B_vals_s[i], koff_0_vals[i])
            koff_fit.append(koff_s)
        
        # catch-slip model
        elif (len(lists) == 7) or (len(lists) == 11):
            koff_cs = catch_slip(f, E_21_vals[i],
                                    k_1rup_vals[i],
                                    f_12_vals[i],
                                    k_2rup_vals[i],
                                    x_B_vals_cs[i])
            koff_fit.append(koff_cs)
    return koff_fit
    
# koff fit vals barely change
koff_fit_vals = koff_fit_func(f_fit_vals)

# force to Pe conversion function
def f_to_Pe_func(f):
    Q = f * np.sqrt((2*L)/a) * ((w*b**2) / (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2))
    tau = (3*mu*Q) / (2*w*b**2)
    shear_rate = tau/mu
    u_f = (y*10**(-6)) * shear_rate * (1 - (5/16)*(a/y)**3) # u_f in units of um/s
    V = u_f * (1-speed_const)
    Pe = (V*alpha)/D # alpha, D have um in units
    return Pe

Pe_fit_vals = f_to_Pe_func(f_fit_vals)

#Probability of interaction function w/ fitted parameters
def P_in_fit(Pe):
    delta = ((alpha**2)*kin_avg)/D
    cap_lambda = lambda_func(Pe)
    P = (cap_lambda*delta)/(1+cap_lambda*delta)
    return P

#k+ fit function
kplus_fit_vals = Nu_func(Pe_fit_vals)*D*m_r*np.pi*Ac_avg*P_in_fit(Pe_fit_vals)
    
#NbNT fit function
def NbNT_fit_func(f, m_l_arr, koff_fit, kplus_fit):
    NbNT_fit = []
    for i in range(len(m_l_arr)):
        NbNT = (kplus_fit * site_densities[i]) / (kplus_fit * site_densities[i] + koff_fit[i])
        NbNT_fit.append(NbNT)
        
    return NbNT_fit

NbNT_fit_vals = NbNT_fit_func(f_fit_vals, site_densities, koff_fit_vals, kplus_fit_vals)

# %% plot: k_off vs f with curve fits

#color cycle scale: setting each set of data + fitted curve to a different color 
colors = iter(cm.rainbow(np.linspace(0, 1, len(site_densities))))
        
# koff vs f plot

for i in range(len(site_densities)):
    if (len(lists) == 4) or (len(lists) == 11):
        plt.figure(0)
        plt.xlabel('Force (pN)')
        plt.ylabel(r'$k_{off} (s^{-1})$')
        plt.title('Bond Dissociation Model: Slip')
        c = next(colors)
        plt.scatter(force_vals_s[i], koff_vals_s[i], color=c,
                  label=r'Data for $m_l = %d$' u' sites/\u03BCm\u00b2'% site_densities[i])
        
        plt.plot(f_fit_vals, koff_fit_vals[i], color=c,
                  label=r'Slip Model for $m_l = %d$' u' sites/\u03BCm\u00b2' % site_densities[i])
            
        # avoid graphing error when only 1 trial is used
        for i in range(len(koff_error_vals)):
            for j in range(len(koff_error_vals[i])):
                if koff_error_vals[i][j] == koff_error_vals[i][j]:
                    plt.errorbar(force_vals_s[i][j], koff_vals_s[i][j],
                                  yerr=koff_error_vals[i][j],
                                  ecolor=c,
                                  capsize=5,
                                  fmt='none')
        
        plt.legend()
        plt.savefig('slip.png', dpi=300)
    
    elif (len(lists) == 7) or (len(lists) == 11):
        plt.figure(1, figsize=(9,7))
        plt.xlabel('Force (pN)')
        plt.ylabel(r'$k_{off} (s^{-1})$')
        plt.title('Bond Dissociation Model: Catch-Slip')
        c = next(colors)
        plt.scatter(force_vals_cs[i], koff_vals_cs[i], color=c,
                    label=r'Data for $m_l = %d$' u' sites/\u03BCm\u00b2'% site_densities[i])
        
        plt.plot(f_fit_vals, koff_fit_vals[i], color=c,
                 label=r'Catch-Slip Model for $m_l = %d$' u' sites/\u03BCm\u00b2' % site_densities[i])
            
        plt.errorbar(force_vals_cs[i], koff_vals_cs[i], yerr=koff_error_vals[i],
                    ecolor=c, capsize=5, fmt='none')
                    
        plt.legend()
        plt.savefig('catch_slip.png', dpi=300)

# %%   plot: Nb/NT vs f and k+ vs f with curve fits
# NbNT vs. F plot
plt.figure(2)
plt.xlabel('Force (pN)')
plt.ylabel(r'Capture Efficiency ($N_b / N_T$)')

colors = iter(cm.rainbow(np.linspace(0, 1, len(site_densities))))

for i in range(len(NbNT_vals)):
    c = next(colors)
    plt.scatter(forces[i], NbNT_vals[i], color=c, 
             label=r'$m_l = %d$' u' sites/\u03BCm\u00b2' % site_densities[i])
    plt.plot(f_fit_vals, NbNT_fit_vals[i], color=c,
             label=r'Best-fit: $m_l = %d$' u' sites/\u03BCm\u00b2' % site_densities[i])
    
    plt.errorbar(forces[i], NbNT_vals[i], yerr=NbNT_error_vals[i],
                 ecolor=c, fmt='none', capsize=5)
    
plt.legend()
plt.savefig('NbNT.png', dpi=300)
    
# k+ vs f plot
plt.figure(3)
plt.xlabel('Force (pN)')
plt.ylabel(r'$k_{+}$' u' (\u03BCm\u00b2/s)')
plt.scatter(forces[0], kplus_vals_avg, color=c, label=r'AVG $m_l$')
plt.plot(f_fit_vals, kplus_fit_vals, color=c, label=r'best-fit: $k_+$')
plt.errorbar(forces[0], kplus_vals_avg, yerr=kplus_errors, 
             ecolor=c, capsize=5, fmt='none')
plt.legend()
plt.savefig('k+.png', dpi=300)

# %% writing RLNEK koff values with trackID to csv file for all inputed conditions
headers = ['Filename', 'koff (1/s)','TRACK_ID']

with open('RLNEK_koff.csv', 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(['Module 2: koff values'])
    wr.writerow(headers)  
    for i in range(len(koff_vals)):
        for j in range(len(koff_vals[i])): 
            for k in range(len(koff_vals[i][j])):
                for m in range(len(koff_vals[i][j][k])):
                    if m == 0:
                        wr.writerow([track_data[i][j][k], koff_vals[i][j][k][m], koff_trackID_vals[i][j][k][m]])
                    else:
                        wr.writerow(['', koff_vals[i][j][k][m], koff_trackID_vals[i][j][k][m]])
    
# %% writing RLNEK summary statistics to csv file
csv_data = []
for i in range(len(site_densities)):
    csv_sublist = [] # [filename, force, koff, ...] for each site density
    for j in range(len(track_data[i])):
        csv_subsub = [track_data[i][j], m_r, site_densities[i], 
                      Q_vals[i][j], forces[i][j], koff_avg_vals[i][j], koff_error_vals[i][j], NbNT_vals[i][j], 
                      NbNT_error_vals[i][j], kplus_vals[i][j], kplus_star_vals[i][j]]
        
        csv_sublist.append(csv_subsub)
    
    csv_data.append(csv_sublist)
    
headers = ['Filename', 'm_r (sites/\u03BCm\u00b2)', 'm_l (sites/\u03BCm\u00b2)', 
           'Q (\u03BCL/hr)', 'Force (pN)', 'koff (1/s)', 'koff SEM', 'Nb/NT', 
           'Nb/NT SEM', 'k+ (\u03BCm\u00b2/s)', 'k+* (\u03BCm\u00b2/s)']

with open('RLNEK.csv', 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(['Module 2'])
    wr.writerow(headers)

    for i in range(len(csv_data)):
        for j in range(len(csv_data[i])):
            wr.writerow(csv_data[i][j])