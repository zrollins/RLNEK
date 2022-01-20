# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:31:15 2020

@author: Allison
"""
import pandas as pd
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt, cm
import os

#user directory/pathname
os.chdir('/User/RLNEK_tests/NT_Optimization_test')
# %% parameters with unit conversions

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

y = a+d


## SET the experimental system parameters to reduce redundancy
#a = 7.6 * 1e6 # cell radius (um)
#d = 0.050 * 1e6 # distance from wall to edge of cell (um)
#L = 36 * 1e3 # bond length (nm)
#w = 800 * 1e6# c hamber width (um)
#b = 37.5 * 1e6# chamber half height (um)
#mu = 6.92e-3 * 1e-13# (dyne/cm^2)


CCD_FPS = int(input('Enter CCD FPS: '))
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

y = a+d

# %% loop that inputs user data
# need data with same flow rate, different site densities
track_data = []
spots_data = []
site_densities = []
shear_rates = []
forces = []
Q_vals = []
while True:
    run = str(input('Enter \"y\" to input data or \"n\" to quit: '))
    
    if run.lower() == 'n':
        break
    
    elif run.lower() == 'y':
        
        # flow rate
        Q_nc = float(input('Enter flow rate (\u03BCL/hr): ')) 
        Q_vals.append(Q_nc)
        Q = Q_nc * (10**27 / 3600) # microliter/h to pm^3/s conversion
        
        m_l = input('For Q = %.2f, enter site densities (sites/\u03BCm\u00b2): ' % Q_nc)
        m_l_str = [val.strip() for val in m_l.split(',')]
        
        m_l_arr = np.zeros(len(m_l_str))
        for i in range(len(m_l_str)):
            m_l_arr[i] = float(m_l_str[i])
            
        site_densities.append(list(m_l_arr))
            
        # tether force
        f = Q * np.sqrt(a/(2*L)) * (1.7005*9*np.pi*mu*a**2 + 0.9440*6*np.pi*mu*a**2) / (w*b**2)
        forces.append(f)
        
        # shear stress
        tau = (3*mu*Q) / (2*w*b**2)
        shear_rate = tau / mu
        shear_rates.append(shear_rate)
        
        # Trackmate files
        track_data_sublist = []
        spots_data_sublist = []
        t_min_sublist = []
        
        for i in range(len(m_l_arr)):
            track_file_name = input('For flow rate = %.2f (\u03BCL/hr) and site density = %f (sites/\u03BCm\u00b2), enter name of "tracks" file(s) from Trackmate: ' % (Q_nc, m_l_arr[i]))
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
                
            spots_file_name = input('For flow rate = %.2f (\u03BCL/hr) and site density = %f (sites/\u03BCm\u00b2), enter name of "spots" file(s) from Trackmate: ' % (Q_nc, m_l_arr[i]))
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
                        
        track_data.append(track_data_sublist)
        spots_data.append(spots_data_sublist)
                 
    else:
        print('Please enter \"y\" or \"n\".')

# %% calculating number of bound cells/spheres, N_b
Nb_vals = []

#3-dimensional nested list format:
##sub = flow rate
##subsub = flow rate, site density
##subsubsub = flow rate, site density, trial
## e.g., Nb_vals[flow rate][site density][trial]

##AVG and SEM is over all inputted trials for a given condition [site density][flow rate]
for m in range(len(track_data)):
    Nb_vals_sub = []
    
    for n in range(len(track_data[m])):
        Nb_vals_subsub = []
        
        # cell velocity filtering
        u_f = y*shear_rates[m]*(1-(5/16)*(a/y)**3) * 1e-6 # convert back to microns
        
        for p in range(len(track_data[m][n])):
            # old Trackmate file format
            # tracks_raw_data = pd.read_csv(track_data[m][n][p])
            # spots_raw_data = pd.read_csv(spots_data[m][n][p])
            
            # new Trackmate format
            tracks_raw_data = pd.read_csv(track_data[m][n][p], header=0,skiprows=range(1,4), encoding= 'unicode_escape')
            spots_raw_data = pd.read_csv(spots_data[m][n][p],  header=0,skiprows=range(1,4), encoding= 'unicode_escape')
            
            filtered_speeds = tracks_raw_data[tracks_raw_data['TRACK_MEAN_SPEED'] < np.absolute(u_f)]
            filtered_tracks_list = list(filtered_speeds['TRACK_ID'])
            
            # only collect track IDs from spots stats file for certain velocities
            better_tracks = []
            
            # obtaining track ID's present in both spots and track stats spreadsheets
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
                            
            # new lists to categorize data after velocity filtering
            particleID_new = []
            trackID_new = []
            x_new = []
            y_new = []
            frame_new = [] 
            
            # adding better_tracks corresponding data to empty lists 
            for i in range(len(better_tracks)):
                for j in range(len(trackID)):
                    if trackID[j] == better_tracks[i]:
                        particleID_new.append(particleID[j])
                        trackID_new.append(trackID[j])
                        x_new.append(x_pos[j])
                        y_new.append(y_pos[j])
                        frame_new.append(frame[j])
            
            # r refers to meeting criteria
            # find which particles meet stopping criteria
            r_pos_x = []
            r_pos_y = []
            r_trackID = []
            r_particleID = []
            r_frame = []
            
            i = 0
            i_max = len(trackID_new)
            j = 0
            
            # calculate displacement of a given cell
            def calc_disp(x0,x,y0,y):
                return np.sqrt((x-x0)**2+(y-y0)**2)
            
            # filter using stopping criteria
            tmin_frames = t_min_input * CCD_FPS
            while i < i_max-1:
                disp1 = calc_disp(x_new[i+1],x_new[j],y_new[i+1],y_new[j])
                if disp1 <= 1:
                    i += 1
                    disp2 = calc_disp(x_new[i],x_new[j],y_new[i],y_new[j])
                    if i-j > tmin_frames:
                        r_particleID.append(particleID_new[i])
                        r_trackID.append(trackID_new[i])
                        r_pos_x.append(x_new[i])
                        r_pos_y.append(y_new[i])
                        r_frame.append(frame_new[i])
                else:
                    i += 1
                    j = i-1
            
            # stopping events time conversion: (# of frames) -> seconds
            # tc = time conversion
            tc_particleID = np.array(r_particleID)
            tc_trackID = np.array(r_trackID)
            tc_frame = np.array(r_frame)
            tc_pos_x =np.array(r_pos_x)
            tc_pos_y =np.array(r_pos_y)
            
            # initial parameters
            t_total = []
            tc_trackID_new = []
            i = 1
            j = 0
            t_tot = 0
            
            # time conversion
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
                    
            #total number of bound cells/spheres
            Nb = len(tc_trackID_unique)
            Nb_vals_subsub.append(Nb)
            
        # AVG number of bound cells/spheres per condition [flow rate][site density] for replicates
        Nb_avg_new = np.mean(Nb_vals_subsub)
        Nb_vals_sub.append(Nb_avg_new)

    Nb_vals.append(Nb_vals_sub)

# N_T optimization over >= 3 site densities
def Nb_func(mrml,NT,AcKa):
    dem = AcKa*mrml
    return NT/(1+(1/dem))

mrml_vals = []
for i in range(len(site_densities)):
    mrml = m_r * np.array(site_densities[i])
    mrml_vals.append(mrml)
    
N_T_vals = []
AcKa_vals = []
for i in range(len(Nb_vals)):
    params, matrix = optimize.curve_fit(Nb_func, mrml_vals[i], Nb_vals[i],
                                        bounds=[[0,0],[np.inf,np.inf]])
    N_T_vals.append(params[0])
    AcKa_vals.append(params[1])
    
for i in range(len(forces)):
    print('')
    print('For flow rate = %.4f \u03BCL/hr & force = %.4f pN: N_T = %.4f.' % (Q_vals[i], forces[i], N_T_vals[i]))

#curve-fit N_T Optimization
mrml_fit_vals = np.linspace(mrml_vals[0][0], mrml_vals[0][-1], 1000)
def Nb_opt_fit(mrml):
    Nb_fit = []
    for i in range(len(Q_vals)):
        Nb_sub = Nb_func(mrml, N_T_vals[i], AcKa_vals[i] )
        Nb_fit.append(Nb_sub)
    return Nb_fit

Nb_fit_vals = Nb_opt_fit(mrml_fit_vals)
# %%plot NT Optimization
colors = iter(cm.rainbow(np.linspace(0, 1, len(Q_vals))))
for i in range(len(Q_vals)):
    c = next(colors)
    plt.scatter(mrml_vals[i], Nb_vals[i], color=c, label='$Q = %.4f$' % (Q_vals[i]))
    plt.plot(mrml_fit_vals, Nb_fit_vals[i], color=c,
             label='$Q = %.4f$' % (Q_vals[i]))
        
    
plt.xlabel(r'$m_{r}m_{l} (sites/\mu m^{4})$')
plt.ylabel(r'Number of Bound Cells/Spheres, $N_b$')
plt.legend()
plt.savefig('NT_Optimization.png', dpi=300, bbox_inches='tight')