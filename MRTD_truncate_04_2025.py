# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 11:39:31 2025

@author: mrc582
"""

# %% imnport the packages
import os
os.chdir('G:\\My Drive\\Projects\\RIR denoising\\rir-truncation\\') # location of src 
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import scipy
# from coherenceFunctions import * # this is the file with functions to calculate short-time coherence and estimate volatility
from coherenceFunctions import *
from onset_detection_functions import *
from utils import*
from scipy.io import savemat

import pyfar as pf
import pyrato as ra


# %% # load the RIRs used for the analysis - here we are using the MRTD dataset
# get the current working directory
current_working_directory = os.getcwd()

# path to MOSAIC RIRs
path_rir_mosaic = "G:\\My Drive\\Projects\\Blind Multi Room\\rirs_med\\rirs_med"#os.path.join(current_working_directory, "rirs_med", "rirs_med")
print(path_rir_mosaic)

# initialize list to store file names
file_list_mosaic = []

# save all the RIR names 
for (root, dirs, files) in os.walk(path_rir_mosaic):
    file_list_mosaic.extend(files)  # Extend instead of overwrite
    
# print(file_list_mosaic)  
num_rirs = np.shape(file_list_mosaic)[0]
print(num_rirs)
# set the path to get single measured RIRs
path_rir_single = "G:\\My Drive\\Projects\\Blind Multi Room\\rirs"#os.path.join(current_working_directory, "rirs")
print(path_rir_single)

# %% declare variables
n_spk = 4 # number of sound sources in the MRTD 
n_ir = 3 # number of repeated RIRs in the MRTD

# for aligning 
upFactor = 10
refID = 0
rirs_aligned = np.zeros([2*48000, n_spk, n_ir+1, 100])
# %% read RIRs and align them
num_rirs = 50
for rir in range(num_rirs):
    # read the mosaic RIRs
    name = file_list_mosaic[rir]
    fs, rir_mosaic = scipy.io.wavfile.read(path_rir_mosaic + '\\' +name)
  
    #read the single RIRs
    rir_single = np.zeros([np.shape(rir_mosaic)[0],n_spk, n_ir])
    for s_rir in range(n_ir):
        _, rir_single[:,:,s_rir] = scipy.io.wavfile.read(path_rir_single + '\\' + name[:-4] + "_ir" + str(s_rir) + ".wav")
        

    rir_len = np.shape(rir_mosaic)[0]
    
    #align the RIRs just in case
    for spk in range(n_spk):
       
        rirs_to_align = np.concatenate((rir_mosaic[:, spk].reshape(rir_len, 1), rir_single[:,spk,:]), axis =1)
        # print(rirs_to_align)
        # first upsample the signal by a factor
        rirs_upsampled = signal.resample(rirs_to_align, (rir_len)*upFactor, axis=0)
        # then time-align sweeps from w.r.t reference one by using cross-correlation
        rirs_aligned_up = time_align_rirs(rirs_upsampled, refID)
        # downsample at the end
        rirs_aligned[:, spk, :, rir]= signal.resample(rirs_aligned_up,rir_len, axis = 0 )
        
# %% declare frequencies for filtering
fcentre  = np.array([250, 500, 1000, 2000, 4000, 8000, 16000])
# fcentre  = np.arange(200, 16000, 200)
fd = 2**(1/2);
fupper = fcentre * fd
flower = fcentre / fd
num_freq = len(fcentre)

num_rirs =50

# %% declare some more
filtered_signal = np.zeros([rir_len, n_ir+1, num_rirs, num_freq]) #for now only for one source, will extend later

# %% and some more
noise_len = int(fs/2)
coh = np.zeros([rir_len,   num_rirs, num_freq])
coh_same = np.zeros([rir_len,  num_rirs, num_freq])
coh_isReg = np.zeros([rir_len,   num_rirs, num_freq])


t_n_is  = np.zeros([ num_rirs, num_freq])
t_s_is  = np.zeros([ num_rirs, num_freq])

var_ref =  np.zeros([  num_rirs, num_freq])
var_both =  np.zeros([  num_rirs, num_freq])
# methods to find onset
M = np.zeros([ num_rirs, num_freq])
M_5 = np.zeros([  num_rirs, num_freq])
D_E = np.zeros([ num_rirs, num_freq])
D_S = np.zeros([  num_rirs, num_freq])
E = np.zeros([ num_rirs, num_freq])

time = np.linspace(0, rir_len / fs, rir_len)
winLen = int(0.025*fs)

wL = 0.01*fs
win = scipy.signal.windows.boxcar(int(wL))
win = win/np.sum(win)

# for lundeby's method 
inter_time  = np.zeros([ num_rirs, num_freq])
rev_time  = np.zeros([ num_rirs, num_freq])
noise_level  = np.zeros([  num_rirs, num_freq])

# %% # Apply the filter
for it in range(num_rirs):
    for it_ in range(n_ir+1):
        for n_freq in range(num_freq):
            filtered_signal[:, it_, it, n_freq] = octave_filter(rirs_aligned[:, 3, it_, it], fs, flower[n_freq], fupper[n_freq]) 
            
# %% truncation point estimation

for n_freq in range(num_freq):  
    print(fcentre[n_freq])
    
   
    for rir in range(num_rirs):
        rirs_analyzed = filtered_signal[:, :,rir,  n_freq]

        for it in range(1):
            # the lundeby method
            pr_rir = pf.Signal(rirs_analyzed[:, it], fs)
            # pr_max = np.abs(pr_rir.time).max()
            # if pr_max > 0:
            #     pr_rir = pr_rir/pr_max
                
            try:
                inter_time[rir,  n_freq], rev_time[rir, n_freq], noise_level[rir, n_freq] = ra.intersection_time_lundeby(pr_rir, fcentre[n_freq])
            except Exception as e:
                print(f"Lundeby method failed for it={it}, n_freq={n_freq}: {e}")
                
                
            
               
            pcc = np.corrcoef(rirs_analyzed, rowvar=False)
            best_m = np.argmax(pcc[1:, 0])
    
            ref_sig = rirs_analyzed[:, it]
            other_sig = rirs_analyzed[:, best_m+1] #np.delete(rirs_analyzed, it, 1)
            
            M[rir,  n_freq] = max_abs_rir(ref_sig)/fs
            M_5[rir,  n_freq]= max_abs_RIR_5ms(M[rir,  n_freq]*fs, fs)/fs
            D_E[rir,  n_freq], e_sig =  mean_over_time(ref_sig, fs)
            D_E[rir,  n_freq] = D_E[rir,  n_freq]/fs
            E[rir,  n_freq] = threshold_E(e_sig) /fs
         
            r_coh = np.convolve(ref_sig*other_sig, win, mode='same')

            coh[:,  rir, n_freq] = r_coh

            r_coh_uni, ind = unimodal_regression(r_coh, time.reshape(-1, 1))
            coh_isReg[:, rir, n_freq] = r_coh_uni

            var_both[ rir, n_freq] = np.var((ref_sig[-noise_len:] +other_sig[-noise_len:])/2) 
            
   
              
            
            vv= var_both[ rir, n_freq]#np.var(ref_sig[-noise_len:])
            _ , t_n_is[ rir, n_freq] = onset_offset_points(r_coh_uni, vv, fs)
            
            var_ref[ rir, n_freq] = np.var(ref_sig[-noise_len:] ) 
            vv= var_ref[ rir, n_freq] 
            t_s_is[ rir,n_freq], _ = onset_offset_points(r_coh_uni, vv, fs)
           

# %% write to matlab for easier plotting
mdic = {"t_s_is": t_s_is, "t_n_is": t_n_is, "var_both": var_both, "M": M, "M_5": M_5, "D_E": D_E, "E": E, "inter_time_mrtd": inter_time, "coh_isReg_mrtd": coh_isReg, "filtered_signal_mrtd": filtered_signal}
savemat(r'G:\My Drive\Projects\RIR denoising\rir-truncation\MRTD_results_2.mat', mdic)