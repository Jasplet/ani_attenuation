#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 09:51:11 2022

@author: ja17375
"""

import numpy as np
import obspy
import numba

@numba.jit()
def apply_tstar_operator(trace, fref, tstar):
    """
    Applies a causal t* operaor (at a reference fruency fref) to the input trace

    Reference: Matheney, M. P. & Nowack, R. L. Geophys J Int 123, 1–15 (1995).
    
    Ported from msac_apply_tstar_operator by J Wookey (2022)
    
    Parameters
    ----------
    trace : TYPE
        DESCRIPTION.
    fref : TYPE
        DESCRIPTION.
    tstar : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    signal = trace.data   
    # Take fft of trace. Supplying a larger n 0-pads trace.
    n = int(nextpow2(trace.stats.npts))
    fd_signal = np.fft.fft(signal, n)
    frequencies = np.fft.fftfreq(n, d=trace.stats.delta)[0:n//2]
    #Negative frequencies are in second half of array     
    
    #get angular frequencies
    ang_freqs = 2*np.pi*frequencies
    ang_fref = 2*np.pi*fref
    # Create causual t* multiplier 
    aw_imag = (-1j*tstar*ang_freqs[1:]/np.pi)*np.log(ang_freqs[1:]/ang_fref)
    aw_real = -1*ang_freqs[1:]*tstar/2
    aw = np.ones(int(n/2,),dtype=np.complex128())
    aw[1:] = np.exp(aw_real + aw_imag)
    # Apply t* operate to frequency domain signal
    attenuated_fd_signal = np.zeros((int(n),), dtype=np.complex128())
    attenuated_fd_signal[0:n//2] = fd_signal[0:n//2]*aw
    # Restore symettry for negative frequencies
    inds = np.arange(n//2-1,-1,-1, dtype=np.int32())
    attenuated_fd_signal[n//2:] = np.conjugate(attenuated_fd_signal[inds])
    # Take inverse fft
    attenuated_signal = np.fft.ifft(attenuated_fd_signal)
    
    trace.date = attenuated_signal[0:trace.stats.npts].real
    return trace

def nextpow2(i):
  n = 2
  while n < i: n = n*2
  return n

if __name__ == '__main__':
    st = obspy.read('/Users/ja17375/Projects/Matisse_Synthetics/ppv1/ideal/Noise20/data/SWAV01.BHE')
    trace = st[0]
    trace.plot()
    apply_tstar_operator(trace, 1, 1)
    trace.plot()