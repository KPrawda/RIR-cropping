# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:09:10 2024

@author: K. Prawda
"""

# coherence functions used in RIR processing
# import relevant packeges
import numpy as np
import scipy
import matplotlib.pyplot as plt
# %%

def slidingCoherence(irRef,ir,winLen):
    # Input:
    # - irRef = reference RIR
    # - ir = RIRs to calculate coherence
    # - winLen = length of the analysis window in samples
    # Output:
    # - r_coh = calculated coherence
    # - e_sig = RIR energy
    # - r_snr = SNR-based expected coherence
    # - e_ref = energy of the reference RIR
    # - e_mean = energy of the averaged RIRs    
    
    if np.shape(np.shape(ir))[0] ==1:
        irLen = np.shape(ir)[0]
        numIR =1
        ir = np.reshape(ir, (irLen, numIR))
    else:
        irLen, numIR = np.shape(ir)
    
    print(numIR)
    # estimate noise levels
    noise = irRef[int(np.round(0.9*irLen)):]
    noiseLevel = np.sqrt(np.median(noise**2,axis = 0))
    # print(noiseLevel)
    
    # window
    win = scipy.signal.windows.boxcar(winLen)
    win = win/np.sum(win)
    
    # calculate coherence    
    e_sig = np.zeros([irLen, numIR])
    e_mean = np.zeros([irLen, numIR])
    r_coh = np.zeros([irLen, numIR])
    e_ref = np.convolve(irRef**2, win, mode='same')
    
        
    for nir in range(numIR):
        r_cov = np.convolve(irRef*ir[:, nir], win, mode='same')
        e_sig[:, nir] = np.convolve(ir[:,nir]**2, win, mode='same')
        r_coh[:, nir] = r_cov/np.sqrt(e_sig[:, nir]*e_ref)
        avg_s = (irRef+ir[:, nir])/2
        e_mean[:,nir] =  np.convolve(np.abs(avg_s)**2, win, mode='same')
    
    # RIR energy
    e_rir = e_sig - noiseLevel**2
    # plt.figure()
    # plt.plot(e_rir)
    # plt.plot(e_sig)
    # plt.show()
    e_rir[e_rir < 0] = 0
    
    # expected coherence
    r_snr = e_rir/e_sig
    
    return r_coh, e_sig, r_snr, e_ref, e_mean, r_cov

# %% coherence model
def coherenceModel(F,T,volatility):
    # model the spectral cross-coherence for RIR analysis
    # input parameters: 
    # - F, vector of frequencies, shape (n,1)
    # - T, vector of time, shape (1,m)
    # - volatility, the value of volatility, so far a scalar
    #  output parameters:
    # - pred_coh, predicted coherence, shape (n,m)
    # - pred_toastd, predicted TOA standard deviation, shape (1,m)
    magicFactor = 20
    pred_coh = np.exp(- magicFactor * np.square(F * np.sqrt(np.transpose(T)) * volatility))
    pred_toastd = T * volatility
    
    return pred_coh, pred_toastd



# %%

def findVolatility(time_vec, meas_coh, mask, snr_coh, fb):
    # fit the volatility to the measured data
    # input parameters: 
    # - time_vec, vector of time for meas_coh, shape (1,m) or (m,)
    # - meas_coh, measured coherence between RIRs, shape (m, n, b)
    # - mask, high energy region, shape (m,n,b)
    # - snr_coh, coherence expected from snr, shape (m,n,b)
    # - fb, center frequencies for bands, shape (b,)
    #  output parameters:
    # - volatility, the value of volatility, shape (b,)
    
    if np.shape(np.shape(meas_coh)) == (1,):
        numT= np.shape(meas_coh)
        numIR =1
        numB =1
    else:
        numT, numIR, numB = np.shape(meas_coh)
    
    # print(numT, numIR, numB)
    
    volatilityMin = -50;
    volatilityMax = -3;
    tolerance = 0.01
    volatility = np.zeros([numIR, numB])
    
    for ir in range(numIR):
        for band in range(numB):
            m = mask[:, band]
            
            T = time_vec[m==1]
            coh = meas_coh[m==1, ir,band]
            F  =fb[band]
            snr  = snr_coh[m==1, ir,band]
            
            def loss_fun(volatility):
                temp,_ = np.abs(coherenceModel(F,T,np.exp(volatility)) * snr - coh)
                return np.sum(temp)
            
            vol_temp= scipy.optimize.fminbound(loss_fun, volatilityMin, volatilityMax, xtol=tolerance)#, full_output='true', disp = 3)
            # print(fval)
            volatility[ir,band] = np.exp(vol_temp)
    
    return volatility
    