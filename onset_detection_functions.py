# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:47:12 2025

@author: mrc582
"""

# onset finding methods
import numpy as np
from scipy import signal
import scipy

# %% maximum of the abs RIR

def max_abs_rir(rir):
    M = np.argmax(np.abs(rir))
    return M

# %% maximum -5 ms

def max_abs_RIR_5ms(M,fs):
    offset = int(np.round(0.005*fs))
    M_5 = M-offset
    return M_5


# %% mean over time D_e

def mean_over_time(rir, fs):
    # window
    winLen = 0.01*fs
    win = scipy.signal.windows.boxcar(int(winLen))
    win = win/np.sum(win)
    
    e_sig = np.convolve(rir**2, win, mode='same')
    
    max_e = np.argmax(e_sig)
    
    ratio = np.zeros([ max_e-1, 1])
    for n in range(max_e-1):
        ratio[n] = e_sig[n+1]/e_sig[n]
        
    D_E = np.argmax(ratio)
    return D_E, e_sig

# %% threshold E

def threshold_E(e_sig):
    K=3
    e_len = len(e_sig)
    E = 0

    for n in range(e_len):
        median_energy = np.median(e_sig[:n])
        if e_sig[n] >= K*median_energy:
            E = n
            break
        
    return E

# %% mean over spectra D_s

def mean_over_spectra(rir, fs):
    # spectrogram parameters
    overlap = 128
    winLen = 2*overlap
    nfft = 2**10
    
    f,t,rir_stft = signal.spectrogram(rir, fs = fs, window = 'hamming', nperseg= winLen, noverlap = overlap, nfft = 2*nfft-1,  return_onesided='true', axis =0)  # perform the STFT
    
    #sum over the freq buns
    E_n = np.sum(rir_stft, axis=0)
    
    max_e = np.argmax(E_n)
    
    ratio = np.zeros([ max_e-1, 1])
    for n in range(max_e-1):
        ratio[n] = E_n[n+1]/E_n[n]
        
    D_S = np.argmax(ratio)
    return D_S

# %% read the RIR and find onset points

def find_onset(rir, fs):
    
    M = max_abs_rir(rir)
    M_5 = max_abs_RIR_5ms(M,fs)
    D_E, e_sig =  mean_over_time(rir, fs)
    E = threshold_E(e_sig)
    # D_S = mean_over_spectra(rir, fs)
    D_S = 0
    return M, M_5, D_E, E, D_S