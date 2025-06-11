# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:56:29 2025

@author: K. Prawda
"""

# %% imports
import os
os.chdir('G:\\My Drive\\Projects\\RIR denoising\\rir-truncation\\') # location of src, if needed 
from coherenceFunctions import *
from utils import*
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.io import savemat


# %% to db
def db(signal, scale):
    # scale is 1 or 2
    signal_db = scale*10*np.log10(np.abs(signal))
    
    return signal_db

# %% read the synthetic RIRs (they have some volatility baked in )
fs, rir_1 =  scipy.io.wavfile.read("IR_Synthetic_1.wav")
fs, rir_2 =  scipy.io.wavfile.read("IR_Synthetic_2.wav")

# normalize
rir_1 = rir_1/np.max(np.abs(rir_1))
rir_2 = rir_2/np.max(np.abs(rir_2))

rir_len = len(rir_1)
time = np.linspace(0, rir_len/fs, rir_len)

## apply decay to RIRs
RT = 1 #s
tau = 6.908/RT
rir_decay_1 = rir_1*np.exp((-tau*time)) #1s decay
rir_decay_2 = rir_2*np.exp(-0.5*tau*time) # 2s decay

rir_decay_1 = rir_decay_1/np.max(np.abs(rir_decay_1))
rir_decay_2 = rir_decay_2/np.max(np.abs(rir_decay_2))

# create random noise sequencyes
rng = np.random.default_rng(42)
noise_1 = rng.standard_normal(rir_len)  
noise_1 = noise_1/np.max(np.abs(noise_1))

rng2 = np.random.default_rng(150)
noise_2 = rng2.standard_normal(rir_len)
noise_2 = noise_2/np.max(np.abs(noise_2))

# declare SNRs
SNR = np.arange(5,35,5)-5
num_cases = len(SNR);
SNR_dec =  [s/20 for s in SNR]
SNR_dec = np.round(SNR_dec, 2)
noise_gain = 10**-SNR_dec

# move the onset by 0.5 s
rir_decay_1 = np.roll(rir_decay_1, int(0.5*fs), axis=0)
rir_decay_2 = np.roll(rir_decay_2, int(0.5*fs), axis=0)

rir_3 = rir_decay_2 - 0.5*rir_decay_1 # fade in 
rir_3 = rir_3/np.max(np.abs(rir_3))
# %% ground-truth onset and truncation points 
# bias = np.arange(1,7, 1)/12 # 30 dB of decay, in 5 dB steps for single slope
bias = np.arange(1,7, 1)/6 # 30 dB of decay, in 5 dB steps for double slope
op_ground_truth = 0.5 # the RIR was rolled by this much
tp_ground_truth=op_ground_truth+bias
# %% correlated, uncorrelated vs anticorrelated

# first declare all the variables
wL = 2**10
win = scipy.signal.windows.boxcar(int(wL))
win = win/np.sum(win)

rir_ref = np.zeros([rir_len, num_cases])
rir_cor = np.zeros([rir_len, num_cases])
rir_unc = np.zeros([rir_len, num_cases])
rir_ant = np.zeros([rir_len, num_cases])
rir_wcr = np.zeros([rir_len, num_cases])

r_cov_cor = np.zeros([rir_len, num_cases])
r_cov_unc = np.zeros([rir_len, num_cases])
r_cov_ant = np.zeros([rir_len, num_cases])
r_cov_wcr = np.zeros([rir_len, num_cases])
r_cov_cor_u = np.zeros([rir_len, num_cases])

v_n_cor = np.zeros([num_cases])
v_n_unc = np.zeros([num_cases])
v_n_ant = np.zeros([num_cases])
v_n_wcr = np.zeros([num_cases])

tp_cor = np.zeros([ num_cases])
tp_unc = np.zeros([ num_cases])
tp_ant = np.zeros([ num_cases])
tp_wcr = np.zeros([ num_cases])
tp_cor_u = np.zeros([ num_cases])

op_cor = np.zeros([ num_cases])
op_unc = np.zeros([ num_cases])
op_ant = np.zeros([ num_cases])
op_wcr = np.zeros([ num_cases])
op_cor_u = np.zeros([ num_cases])


cr = 0.2 # how correlated the noise sequences need to be

noise_3  = cr * noise_1 + np.sqrt(1 - cr**2) * noise_2 # create correlated noise 

# %% the main loop 
norm_factor = np.zeros([ num_cases])
for n in range(num_cases):
    
    # uncomment for the single slope and no fade-in
    # rir_ref[:, n] = rir_decay_1 + noise_1 * noise_gain[n] # correlated noise
    # rir_unc[:, n] = rir_decay_1 + noise_2 * noise_gain[n] # uncorrelated noise
    # rir_ant[:, n] = rir_decay_1 - noise_1 * noise_gain[n] # anticorrelated noise
    # rir_wcr[:, n] = rir_decay_1 + noise_3 * noise_gain[n] # weekly correlated noise
    
    # uncomment for the fade-in
    rir_ref[:, n] = rir_3 + noise_1 * noise_gain[n] # correlated noise
    rir_unc[:, n] = rir_3 + noise_2 * noise_gain[n] # uncorrelated noise
    rir_ant[:, n] = rir_3 - noise_1 * noise_gain[n] # anticorrelated noise
    rir_wcr[:, n] = rir_3 + noise_3 * noise_gain[n] # weekly correlated noise
    
      
    _, _,_,_,_,r_cc = slidingCoherence(rir_ref[:, n],rir_ref[:, n],2*wL)
    _, _,_,_,_,r_c = slidingCoherence(rir_ref[:, n],rir_unc[:, n],2*wL)
    _, _,_,_,_,r_ca = slidingCoherence(rir_ref[:, n],rir_ant[:, n],2*wL)
    _, _,_,_,_,r_cw = slidingCoherence(rir_ref[:, n],rir_wcr[:, n],2*wL)

    r_ca=r_ca-np.min(r_ca)
    norm_factor[n] =np.max(np.abs(r_cc)) # normalization for noise variance
    
    #normalize the covariance curves
    r_cc = r_cc/np.max(np.abs(r_cc))
    r_c = r_c/np.max(np.abs(r_c))
    r_ca = r_ca/np.max(np.abs(r_ca))
    r_cw = r_cw/np.max(np.abs(r_cw))
    
    
    r_cov_cor[:, n], _ = unimodal_regression(np.abs(r_cc), time)
    r_cov_unc[:, n], _ = unimodal_regression(np.abs(r_c), time)
    r_cov_ant[:, n], _ = unimodal_regression(np.abs(r_ca), time)
    r_cov_wcr[:, n], _ = unimodal_regression(np.abs(r_cw), time)

   
    # noise variance
    v_n_cor[ n] = np.var(noise_1 * noise_gain[n])/norm_factor[n]
    v_n_unc[ n] = np.var((noise_1 * noise_gain[n]+ noise_2 * noise_gain[n])/2)
    v_n_ant[ n] = np.var((noise_1 * noise_gain[n]- noise_1 * noise_gain[n])/2)
    v_n_wcr[ n] = np.var((noise_1 * noise_gain[n]+ noise_3 * noise_gain[n])/2)
    

    # estimate onset and trucation points
    op_cor[ n], tp_cor[ n] = onset_offset_points(r_cov_cor[:, n], v_n_cor[ n], fs)
    op_unc[ n], tp_unc[ n] = onset_offset_points(r_cov_unc[:, n], v_n_cor[ n], fs)
    op_ant[ n], tp_ant[ n] = onset_offset_points(r_cov_ant[:, n], v_n_cor[ n], fs)
    op_wcr[ n], tp_wcr[ n] =onset_offset_points(r_cov_wcr[:, n], v_n_cor[ n], fs)
    
    # another option - add uncorrelated noise to fully correlated noise
    _, _,_,_,_,r_cu = slidingCoherence(rir_ref[:, n],rir_ref[:, n]+ noise_2 * noise_gain[n],2*wL)
    r_cu = r_cu/np.max(np.abs(r_cu))
    r_cov_cor_u[:, n], _ = unimodal_regression(np.abs(r_cu), time)
    op_cor_u[ n], tp_cor_u[ n] = onset_offset_points(r_cov_cor_u[:, n], v_n_cor[ n], fs)


# %% noise variance for different noise portion lengths
v_noise_len = np.zeros([num_cases, num_cases])
tp_unc_len = np.zeros([num_cases, num_cases])
op_unc_len = np.zeros([num_cases, num_cases])
lens = [ 0.1*fs, 0.2*fs, 0.4*fs, 0.5*fs, fs, 0.7*fs]
lens = [int(L) for L in lens]
for n in range(num_cases):
    for n2 in range(num_cases):
        v_noise_len[n,n2] = np.var(noise_1[:lens[n]] * noise_gain[n2])/norm_factor[n2]
        op_unc_len[ n,n2], tp_unc_len[ n, n2] = onset_offset_points(r_cov_unc[:, n2], v_noise_len[n, n2], fs)
#%% 

n = 0
f, axs = plt.subplots(2,3, sharey=True, sharex=True)
f.set_figheight(6)
f.set_figwidth(16)
for ax, n in zip(axs.flat, range( num_cases)):
    ax.plot(time, r_cov_cor[:, n], 'k', label = 'RIRs + corr. noise')
    ax.plot(time, r_cov_unc[:, n], 'r', label = 'RIRs + uncorr. noise')
    ax.plot(time, r_cov_ant[:, n], 'c', label = 'RIRs + anticorr. noise')
    ax.plot(time, r_cov_wcr[:, n], 'm', label = 'RIRs + weakly corr. noise')
    ax.plot(time, r_cov_cor_u[:, n], 'g--', label = 'RIRs + cor + unc noise')
    
    
    
    ax.hlines( v_n_cor[ n], 0,2, 'g', label ='Noise variance')
    # ax.vlines(t_points[n], -1, 1, color = 'gold', label = '0.25 coherence')
    ax.vlines(tp_ground_truth[n], -1, 1, color = 'gold', label = 'GT')
    # plt.plot(time, r_n_cor[:, n], 'k--', label = 'correlated noise')
    # # plt.plot(time, 2*r_n_cor[:, n])
    # plt.plot(time, r_n_unc[:, n], 'r--', label = 'uncorrelated noise')
    # plt.plot(time, r_n_ant[:, n], 'c--', label = 'anticorrelated noise')
    # plt.plot(time, r_n_wcr[:, n], 'm--', label = 'weakly correlated noise')
    
    # ax.set_ylim([10**-7, 2*10**-2])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Covariance (product)')
 
    ax.set_yscale('log')
    
ax.legend(ncol=1)
# %% subtract noise covariance

n = 4
plt.figure()
plt.plot(time, r_cov_cor[:, n]-2*r_n_cor[:, n], 'k', label = 'RIRs - correlated noise')
plt.plot(time, r_cov_unc[:, n]-r_n_cor[:, n], 'r', label = 'RIRs - uncorrelated noise')
plt.plot(time, r_cov_ant[:, n]-r_n_ant[:, n], 'c', label = 'RIRs - anticorrelated noise')
plt.plot(time, r_cov_wcr[:, n]-r_n_cor[:, n], 'm', label = 'RIRs - weakly correlated noise')
plt.vlines(t_points[n], -1, 1, 'k')


ax = plt.gca()
# ax.set_ylim([-0.00025, 0.00025])
ax.set_xlabel('Time [s]')
ax.set_ylabel('Covariance (product)')
ax.legend(ncol=2)
# %% save to mat file for easier plotting later on
mdic = {"tp_ground_truth": tp_ground_truth,"op_ground_truth": op_ground_truth,
        "r_cov_cor": r_cov_cor,"r_cov_unc":r_cov_unc,"r_cov_ant":r_cov_ant,
        "r_cov_wcr": r_cov_wcr, "r_cov_cor_u": r_cov_cor_u,
        "v_n_cor": v_n_cor, "v_n_unc": v_n_unc,"v_n_ant":v_n_ant,  
        "v_n_wcr": v_n_wcr,"op_cor": op_cor, "tp_cor": tp_cor, "op_unc": op_unc, "tp_unc": tp_unc, 
        "op_wcr": op_wcr, "tp_wcr": tp_wcr, "op_ant": op_ant, "tp_ant": tp_ant,
        "op_cor_u": op_cor_u, "tp_cor_u": tp_cor_u, 
        "v_noise_len":v_noise_len, "op_unc_len": op_unc_len, "tp_unc_len": tp_unc_len,"r_c": r_c, "fs":fs}
savemat("simulation_results_fade_in.mat", mdic)
# savemat("simulation_results_no_fade.mat", mdic)