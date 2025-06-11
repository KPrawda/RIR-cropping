# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:39:48 2025

@author: mrc582
"""

# %% imnport the packages
import sys
import os
os.chdir('G:\\My Drive\\Projects\\RIR denoising\\rir-truncation\\') # location of src 
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import scipy
from scipy.signal import butter, filtfilt

from sklearn.isotonic import IsotonicRegression
# from coherenceFunctions import * # this is the file with functions to calculate short-time coherence and estimate volatility
from coherenceFunctions import *
from scipy.io import savemat
import pyfar as pf
import pyrato as ra
from scipy.optimize import isotonic_regression

# %% 
def time_align_rirs(upsampled_rir, referenceID):
    # time-align upsampled sweeps
    num_repeat = upsampled_rir.shape[-1]

   
    for rep in range(referenceID+1, num_repeat):
        # calculate the cross-correlation of upsampled sweeps
        corr = signal.correlate(upsampled_rir[ :,  referenceID], upsampled_rir[ :,  rep])
        # get the time lag as well
        lags = signal.correlation_lags(len(upsampled_rir[:,  referenceID]), len(upsampled_rir[:,  rep]))
        # find where the zero lag is - the lag is symmetric on both sides
        lag_ind = np.nonzero(lags==0)
        zero_lag = lag_ind[0][0] # the position of zero lag
        
        # find the lag of the max correlation - this will define the amount of samples we need to shift
        max_lag = np.argmax(corr)
        sign = np.sign(max_lag - zero_lag) # sign will tell us whether to shift to left or right
        shift = np.abs(max_lag- zero_lag)+1 # needs to be shifted by 1 because of indexing through zero
        
        # shift the IR
        shifted = np.roll(upsampled_rir[ :,  rep], sign*shift)
        upsampled_rir[ :,  rep]= shifted

    return upsampled_rir

# %% 
# 2. Apply an Octave Filter (e.g., filter between 440 Hz and 880 Hz)
def octave_filter(signal, sr, low_freq, high_freq):
    nyquist = 0.5 * sr
    low = low_freq / nyquist
    high = high_freq / nyquist
    b, a = butter(4, [low, high], btype='band')  # 4th order Butterworth bandpass filter
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal

# %%
def unimodal_regression(sig, time):
    """Separate the signal into two regions of increasing and decreasing energy, applyisotonic regression to each of them"""
    """Separate the signal into two regions of increasing and decreasing energy, applyisotonic regression to each of them"""
    """ input:
            sig is the signal vector
            time is the time vector
        output:
            sig_unimodal signal smoothed using unimodal regression
            ind index of the max values in signal -> the increasing regression turns to decreasing
       """
    sig = np.squeeze(sig)
    # we need to separate the increasing and descreasing parts of the signal
    # first, find a maximum point of the signal - we assume there are clearly separable increasing and decreasing sections
    ind = int(np.argmax(sig))
    # ind = ind.astype(int)

    # decreasing regression
    sig_ = sig[ind:]
    time_  = time[ind:]
    # isReg = IsotonicRegression(y_min = 0, y_max = 1, increasing = False, out_of_bounds="clip")
    # sig_isReg = isReg.fit_transform(time_, sig_)
    
    result = isotonic_regression(np.squeeze(sig_), increasing = False)
    sig_isReg_dec = result.x

    if ind > 0:
        #increasing regression
        # isReg = IsotonicRegression(y_min = 0, y_max = 1, increasing = True, out_of_bounds="clip")
        sig_ = sig[:ind]
        time_  = time[:ind]
        # sig_isReg = isReg.fit_transform(time_, sig_)
        
        result = isotonic_regression(sig_, increasing = True)
        sig_isReg_inc = result.x
        sig_unimodal = np.concatenate(( sig_isReg_inc, sig_isReg_dec))
    else: 
        sig_unimodal = sig_isReg_dec

    return sig_unimodal, ind
 
# %% load the RIRs used for the analysis - here we are using the Arni FA dataset
# get the current working directory
current_working_directory = os.getcwd()

# path to FA RIRs
path_rirs = 'G:\My Drive\Projects\Forum Acusticum 2023\RIRs' #os.path.join(current_working_directory, "RIRs_FA")
print(path_rirs)

# initialize list to store file names
file_list = []

# save all the RIR names 
for (root, dirs, files) in os.walk(path_rirs):
    file_list.extend(files)  # Extend instead of overwrite
    
# print(file_list)  
num_rirs = 600# np.shape(file_list)[0]
# %% read and align rirs

upFactor = 10   # upsampling factor for aligning
refID = 0       # which rir should be used as a reference for time aligning

rir_single = np.zeros([110251, num_rirs])
for ind in range( num_rirs):#num_rirs):
    # read the RIRs
    name = file_list[ind]
    fs, temp= scipy.io.wavfile.read(path_rirs + '\\' +name)
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

inter_time  = np.zeros([num_rirs, num_freq])
rev_time  = np.zeros([num_rirs, num_freq])
noise_level  = np.zeros([num_rirs, num_freq])

# %% Apply the filter
for it in range(num_rirs):
    for n_freq in range(num_freq):
        filtered_signal[:, it, n_freq] = octave_filter(rirs_aligned[:int(rir_len/2), it], fs, flower[n_freq], fupper[n_freq])
 
# %% run the main loop
for n_freq in range(num_freq):  
    print(fcentre[n_freq])
    temp_sig = filtered_signal[:, :, n_freq]
    # edc = ra.energy_decay_curve_truncation(temp_sig)
    for it in range( num_rirs-1):#:num_rirs-1):
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
        if it in {177, 178, 187,188}:
            continue
        vv= np.var(filtered_signal[-nL:, it, n_freq])
        diff = cov_isReg[:, it,n_freq]- vv#np.mean(n_cov[:, (it, it+1)], axis=1)
        sign_changes = np.where(np.diff(np.sign(diff)))[0]
        c_max = np.argmax(cov[:, it,n_freq])
        sign_changes = sign_changes[sign_changes > c_max]
        t_n[it, n_freq] = sign_changes[0]/fs              

# %% write to matlab for easier plotting
mdic = {"t_n": t_n, "inter_time": inter_time, "filtered_signal": filtered_signal, "cov_isReg": cov_isReg, "fs":fs}
savemat("Arni_results.mat", mdic)