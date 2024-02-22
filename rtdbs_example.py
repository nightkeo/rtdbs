# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 19:56:07 2024

@author: user
"""

from rtdbs.TS_read import TS_read
from rtdbs.rt_main import rt_main

import numpy as np
import os
from shutil import rmtree

#paths
TS_path='.//input_data//TS'
eqdsk_path='.//input_data//geqdsk'
output_path='.//output_data'

#input data for wave
###############################################################################

alpha=1.217 # for a0=a0*alpha - extend magnetic surfaces at limiter

r_start_vect=[1.1,0.0,0.0] #wave start point

n_start_vect=[
    -0.99452189536827329,
    +0.040842524411819833,
    -0.096218957761800145
    ]

R_start=0.62 #start wave from Z=0 at limiter

# frequency_list=[60] #GHz
frequency_list=[20,29,50,55,60,65,70,75] #GHz

moda=0
InversT=0

###############################################################################

shot_number=40458

filename_ts_ne=TS_path+'//'+str(shot_number)+'//'+str(shot_number)+'_n(R)_old.csv'
data_dens=TS_read(filename_ts_ne)

num_ts_list=[23]

for num_ts in num_ts_list:
    for frequency in frequency_list:
        #time from TS
        time_ts=data_dens['time'][num_ts] #s
        r_ts=data_dens['r'][num_ts] #m
        z_ts=np.zeros(len(r_ts)) #m
        ne_ts=data_dens['array'][num_ts]/1e13 #10^19 m^(-3)
        dne_ts=data_dens['darray'][num_ts]/1e13 #10^19 m^(-3)
        
        #create folder for output data
        output_path_all=output_path+'//'+str(shot_number)+'_'+str(int(time_ts*1e6))+'_'+str(frequency)
        if os.path.isdir(output_path_all):
            rmtree(output_path_all, ignore_errors=True)
        os.mkdir(output_path_all)
        
        #time for equ
        # time_eqdsk=int(time_ts*1e3) #ms
        # filename_equ=eqdsk_path+'//'+str(shot_number)+'//'+'g0'+str(shot_number)+'_00'+str(time_eqdsk)+'000'
        time_eqdsk=int(time_ts*1e3) #ms
        filename_equ=eqdsk_path+'//'+str(shot_number)+'//'+'g0'+str(shot_number)+'.00'+str(time_eqdsk) 
        
        #START!
        rt_main(
                alpha,r_start_vect,n_start_vect,R_start,frequency,moda,InversT,
                filename_equ,output_path_all,
                r_ts,z_ts,ne_ts,dne_ts
                )






























































