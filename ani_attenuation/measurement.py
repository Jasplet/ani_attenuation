#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 16:33:37 2022

@author: ja17375
"""

import numpy as np
import numba
import scipy.signal
import obspy

from waveform_tools import attenuate_traces, rotate_traces

def measure_dtstar(files, snr_max, nfast=181, ndts=81):
    
    #loop over array of sac file names and read
    difr_stack = np.zeros((nfast, ndts))
    for file in files:
        waveform_data = obspy.read(file)
        waveforms = window_traces(waveform_data)
        difr_grid = dtstar_gridsearch(waveforms, nfast, ndts)
        difr_stack += difr_grid
    # normalise difr_grid by number of files
    nfiles = len(files)
    difr_stack = difr_stack/nfiles
    # find min value
    grid_shape = (nfast, ndts)
    idx_min = np.unravel_index(np.argmin(difr_stack),shape=grid_shape)

#@numba.jit(parallel=True)
def dtstar_gridsearch(waveforms, nfast, ndts, fref=1):
    

    fast_directions = np.linspace(-90,90,nfast)
    dtstars = np.linspace(0,4.0,ndts)
    trN = waveforms.select(channel='BHN')[0]
    trE = waveforms.select(channel='BHE')[0]
    difrs = np.zeros((nfast, ndts))
    for i in numba.prange(0,nfast):
        trF, trS = rotate_traces(trN, trE ,fast_directions[i])
        # As we only attenaute the fast trace we can caluclate the inst. freq.
        # for the slow trace now        
        inst_freq_trS = measure_inst_freq(trS)
        for j in numba.prange(0, ndts):
            attenuate_traces(trF, fref, dtstars[j])
            inst_freq_trF = measure_inst_freq(trF)
            difrs[i,j] = np.abs(inst_freq_trF - inst_freq_trS)
    
    return difrs
            
#@numba.jit()
def measure_inst_freq(trace):
    '''
    Instantaneous frequency analysis 
    '''
   
    hilbert_trace = scipy.signal.hilbert(trace.data)
    real_trace = np.real(hilbert_trace)
    imag_trace = np.imag(hilbert_trace)
    # Calculate instantaneous amplitude, phase, frequency followign scipy docs
    inst_amplitude = np.abs(hilbert_trace)
    #inst_phase = np.atan(imag_trace/real_trace)
    delta = trace.stats.delta
    imag_gradient = np.gradient(imag_trace, delta)
    real_gradient = np.gradient(real_trace, delta)
    e2 = 1e-5*np.square(inst_amplitude.max()) # damping factor
    inst_frequency = (1/(2.0*np.pi))*((real_trace*imag_gradient - imag_trace*real_gradient)/(np.square(inst_amplitude) + e2)) 
    # Now calculate weighted instantaneous frequency 
    weighted_inst_f = np.sum(inst_frequency*np.square(inst_amplitude)) / np.sum(np.square(inst_amplitude))
                                      
    return weighted_inst_f

#@numba.jit()
def window_traces(st):
    
    st_func = st.copy()
    times = st_func[0].times() + st_func[0].stats.sac.b
    wbeg = st_func[0].stats.sac.a
    wend = st_func[0].stats.sac.f
    window = (times >= wbeg) & (times <= wend)
    for trace in st_func:
        trace.data = trace.data[window]
        
    return st_func

if __name__ == '__main__':
    from synthetics import gen_synthetic_split
    
    st = gen_synthetic_split(30, 1, dtstar=1)
    trF, trS = rotate_traces(st[0], st[1] ,30)
    trSA = attenuate_traces(trS, 1, 1)
    st = window_traces(st)
    difr= dtstar_gridsearch(st, 100, 100)