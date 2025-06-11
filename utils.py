# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 13:57:07 2025

@author: mrc582
"""
# %% imnport the packages
import numpy as np
from scipy import signal
from sklearn.isotonic import IsotonicRegression
from scipy.optimize import isotonic_regression
from scipy.signal import butter, filtfilt
# %% function for time-aligning the RIRs
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

# %% 2. Apply an Octave Filter (e.g., filter between 440 Hz and 880 Hz)
def octave_filter(signal, sr, low_freq, high_freq):
    nyquist = 0.5 * sr
    low = low_freq / nyquist
    high = high_freq / nyquist
    b, a = butter(4, [low, high], btype='band')  # 4th order Butterworth bandpass filter
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal

# %% function for unimodal regression
def unimodal_regression(sig, time):
    """Separate the signal into two regions of increasing and decreasing energy, applyisotonic regression to each of them"""
    """ input:
            sig is the signal vector
            time is the time vector
        output:
            sig_unimodal signal smoothed using unimodal regression
            ind index of the max values in signal -> the increasing regression turns to decreasing
       """
    

    # we need to separate the increasing and descreasing parts of the signal
    # first, find a maximum point of the signal - we assume there are clearly separable increasing and decreasing sections
    ind = int(np.argmax(sig))
    # ind = ind.astype(int)

    # decreasing regression
    sig_ = sig[ind:]
   
    result = isotonic_regression(sig_, increasing = False)
    sig_isReg_dec = result.x

    if ind > 0:
        #increasing regression
        sig_ = sig[:ind]
        
        result = isotonic_regression(sig_, increasing = True)
        sig_isReg_inc = result.x
        sig_unimodal = np.concatenate(( sig_isReg_inc, sig_isReg_dec))
    else: 
        sig_unimodal = sig_isReg_dec

    return sig_unimodal, ind

# %% determine onset and offset points
def onset_offset_points(r_cov, noise_v, fs):
    """Get the values of the onset and offset points w.r.t. noise covariance threshold."""
   
    # get the difference between the covariance and noise variance             
    diff = r_cov - noise_v
    # check for the sign changes
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    # see where is the maximal point for covariance values
    c_max = np.argmax(r_cov)
    # sign changes AFTER the max value indicate the truncation point
    sign_changes_ = sign_changes[sign_changes > c_max]
    trunc_point = sign_changes_[0]/fs  
    
    # sign changes BEFORE the max value indicate the onset point    
    sign_changes_ = sign_changes[sign_changes <= c_max]
    if len(sign_changes_) == 0:
        onset_point=0
    else:
        onset_point=sign_changes_[-1]/fs
    
    return onset_point, trunc_point