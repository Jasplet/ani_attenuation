#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 16:33:37 2022

@author: ja17375
"""

import numpy as np
import numba
import scipy.signal

@numba.jit(parallel=True)
def dtstar_gridsearch(trN, trE):
    
    fast_directions = np.arange(-90,91,1)
    dtstars = np.arange(0,4.05,0.05)
    
    for i in numba.prange(0,len(fast_directions)):
        trF, trS = rotate(trN, trE, fast_directions[i])
        # As we only attenaut the fast trace we can caluclate the inst. freq. for the slow trace now        
        inst_freq trS = measure_inst_freq(trS)
        for dtstars in numba.prange(0, len(dtstars)):
            inst_freq_trF = measure_inst_freq(trF)
            difr = np.abs(inst_freq_trF - inst_freq_trS)
            
    
@numba.jit()
def measure_inst_freq(trace):
    '''
    Instantaneous frequency analysis 
    '''
   
    hilbert_trace = scipy.signal.hilbert(trace.data)
    real_trace = np.real(hilbert_trace)
    imag_trace = np.imag(hilbert_trace)
    # Calculate instantaneous amplitude, phase, frequency followign scipy docs
    inst_amplitude = np.abs(hilbert_trace)
    inst_phase = np.atan(imag_trace/real_trace)
    
    delta = trace.stats.delta
    imag_gradient = np.gradient(imag_trace, delta)
    real_gradient = np.gradient(real_trace, delta)
    e2 = 1e-5*np.square(inst_amplitude.max()) # damping factor
    inst_frequency = (1/(2.0*np.pi))*((real_trace*imag_gradient - imag_trace*real_gradient)/(np.square(inst_amplitude) + e2)) 
    # Now calculate weighted instantaneous frequency 
    weighted_inst_f = np.sum(inst_frequency*np.square(inst_amplitude)) / np.sum(np.square(inst_amplitude))
                                      
    return weighted_inst_f

def window_traces(st, wbeg, wend):
    
    pass
