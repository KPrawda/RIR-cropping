# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:39:48 2025

@author: mrc582
"""

# %% imnport the packages
import sys
import os

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import scipy
from scipy.signal import butter, filtfilt

# from coherenceFunctions import * # this is the file with functions to calculate short-time coherence and estimate volatility
from coherenceFunctions import *
from utils import*
from scipy.io import savemat
import pyfar as pf
import pyrato as ra


# %% load the RIRs used for the analysis - here we are using a subset of the Arni dataset
# get the current working directory
current_working_directory = os.getcwd()

# path to Arni RIRs
src_path = os.path.abspath(os.path.join(current_working_directory, '..', 'RIR-cropping', 'RIR subset', 'Arni'))
sys.path.insert(0, src_path)

# initialize list to store file names
file_list = []

# save all the RIR names 
for (root, dirs, files) in os.walk(src_path):
    file_list.extend(files)  # Extend instead of overwrite
    
# print(file_list)  
num_rirs =  np.shape(file_list)[0] # use only 10 RIRs for brevity of calculations
# %% read and align rirs

upFactor = 10   # upsampling factor for aligning
refID = 0       # which rir should be used as a reference for time aligning

rir_single = np.zeros([110251, num_rirs])
for ind in range( num_rirs):
    # read the RIRs
    name = file_list[ind]
    fs, temp= scipy.io.wavfile.read(src_path + '\\' +name)
    rir_single[:np.shape(temp)[0], ind]=temp
  
   
rir_len = np.shape(rir_single)[0]
#align the RIRs just in case they appear at different times in the audio file. 

# first upsample the signal by a factor
rirs_upsampled = signal.resample(rir_single, (rir_len)*upFactor, axis=0)
# then time-align sweeps from w.r.t reference one by using cross-correlation
rirs_aligned_up = time_align_rirs(rirs_upsampled, refID)
# downsample at the end
rirs_aligned= signal.resample(rirs_aligned_up,rir_len, axis = 0 )

# %% declare frequencies for octave filtering
fcentre  = np.array([250, 500, 1000, 2000, 4000, 8000])
fd = 2**(1/2);
fupper = fcentre * fd
flower = fcentre / fd
num_freq = len(fcentre)

winLen = int(0.025*fs)*np.array([6,5,4,3, 2,1])
time = np.linspace(0, rir_len / (2*fs), int(rir_len/2))

# %% initialize filtered signal variable
filtered_signal = np.zeros([int(rir_len/2), num_rirs, num_freq])

# %% initalize other variables
cov = np.zeros([int(rir_len/2), num_rirs-1, num_freq])

inter_time  = np.zeros([num_rirs-1, num_freq])
rev_time  = np.zeros([num_rirs-1, num_freq])
noise_level  = np.zeros([num_rirs-1, num_freq])

# %% Apply the filter
for it in range(num_rirs):
    for n_freq in range(num_freq):
        filtered_signal[:, it, n_freq] = octave_filter(rirs_aligned[:int(rir_len/2), it], fs, flower[n_freq], fupper[n_freq])
 
# %% run the main loop
for n_freq in range(num_freq):  
    print(fcentre[n_freq])
    temp_sig = filtered_signal[:, :, n_freq]
    
    for it in range( num_rirs-1):
        if it in {177, 178, 187,188}:
            continue
        
        # the lundeby method, using pyfar
        pr_rir = pf.Signal(temp_sig[:, it], fs)
        pr_rir = pr_rir/np.abs(pr_rir.time).max()
        try:
            inter_time[it, n_freq], rev_time[it, n_freq], noise_level[it, n_freq] = ra.intersection_time_lundeby(pr_rir, fcentre[n_freq])
        except Exception as e:
            print(f"Lundeby method failed for it={it}, n_freq={n_freq}: {e}")
       
        # choose a pair of rirs for analysis
        ref_sig = temp_sig[:, it] 
        other_sig =temp_sig[:, it+1]
        # calculate the covariance
        _, _,_,_,_,r_cov = slidingCoherence(ref_sig,other_sig,winLen[n_freq])
        # and store it
        cov[:,  it, n_freq] = np.squeeze(r_cov)
    
# %% apply unimodal regression
cov_isReg = np.zeros([int(rir_len/2), num_rirs-1, num_freq])
for n_freq in range(num_freq): 
    for it in range( num_rirs-1):
        cov_isReg[:, it, n_freq], ind = unimodal_regression(cov[:,  it, n_freq], time.reshape(-1, 1))
# %% define noise variance and calculate the truncation time
nL =8820
t_n = np.zeros([num_rirs-1,num_freq])
for n_freq in range(num_freq): 
    for it in range(num_rirs-1):
        if it in {177, 178, 187,188}: # those are corrupted RIRs
            continue
        vv= np.var(filtered_signal[-nL:, it, n_freq])
        _, t_n[it, n_freq] = onset_offset_points(cov_isReg[:, it, n_freq], vv, fs) # only truncation for Arni as all the RIRs are already starting at t=0
            

# %% write to matlab for easier plotting
# mdic = {"t_n": t_n, "inter_time": inter_time, "filtered_signal": filtered_signal, "cov_isReg": cov_isReg, "fs":fs}
# savemat("Arni_results_compact.mat", mdic)