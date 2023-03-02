#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 09:51:11 2022

@author: ja17375
"""

import numpy as np
import obspy

from scipy.interpolate import PchipInterpolator
from scipy.fft import fft, ifft, fftfreq, next_fast_len

def rotate_traces(y_trace, x_trace, theta):
    '''
    Rotate two obspy traces to the desired rotation angle theta. 
    
    Adapted port of msac_rotate 
    Parameters
    ----------
    y_trace : obspy Trace
        y-direction trace (e.g. North component)
    x_trace : obspy Trace
        x-direction trace (e.g. East component)
    theta : float
        Rotation angle in degrees
        
    Returns
    -------
    trace1p : obspy Trace
        rotated trace1
    trace2p : obspy Trace
        rotated trace2
    '''
    ydata = y_trace.data.copy()
    xdata = x_trace.data.copy()
    y_trace1p = y_trace.copy()
    x_trace2p = x_trace.copy()
    
    ydatap, xdatap = rotate(ydata, xdata, theta)
    y_trace1p.data = ydatap
    x_trace2p.data = xdatap
    y_trace1p.stats.sac['cmpaz'] = y_trace.stats.sac['cmpaz'] + theta 
    x_trace2p.stats.sac['cmpaz'] = x_trace.stats.sac['cmpaz'] + theta
    return y_trace1p, x_trace2p
    
def rotate(y, x, theta):
    """
    Rotate 2 orthogonal traces

    Parameters
    ----------
    y : Numpy array
        array holding y component data (e.g. N component in geog. frame)
    x : Numpy array
        array holding x component data (e.g. E component in geog. frame)

    Returns
    -------
    None.

    """
    if len(y) != len(x):
        raise ValueError('Traces must be same length!')
    if (len(y.shape) > 1 ) or (len(x.shape) > 1):
        raise ValueError('Traces must be 1-D arrays!')
    
    trs = np.vstack((x, y))
    # form rotation matrix
    radtheta = np.deg2rad(theta)
    c = np.cos(radtheta)
    s = np.sin(radtheta)
    rotmat = np.array([
        [ np.cos(radtheta), -np.sin(radtheta)],
        [np.sin(radtheta), np.cos(radtheta)]
        ])
    xp, yp = np.dot(rotmat, trs)
    
    return yp, xp
    
def time_base(delta, nsamps):
    '''
    Creates dummy time base centered on 0. Used to generate wavelet 
    
    use floor div to handle cases where w have an odd number of samples
    Parameters
    ----------
    delta : TYPE
        DESCRIPTION.
    nsamps : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    time =  delta*(np.arange(0,nsamps) - nsamps //2)
    return time 


def apply_tshift(trace, tlag, delta, starttime):
    """
    Ported version of msac_tshft by J Wookey (2011)
    """
    npts_shift = tlag/ delta
    npts = len(trace)
    if npts_shift > npts:
        raise ValueError('Requested shift is larger than trace')
    
    # make zero-padded trace to facilitate tshift
    tr_pad = np.zeros((3*npts,))
    tr_pad[npts:npts*2] = trace 
    time_interp = starttime + (delta * np.arange(0,npts)) - tlag
    time_pad = starttime + (delta * np.arange(-npts,2*npts))
    pchip_interp = PchipInterpolator(time_pad, tr_pad)

    tr_shifted = pchip_interp(time_interp)
    return tr_shifted

def randn_noise(n, amp):
    """
    Generate gaussian random noise trace scald by a given amplitude (should be max amp * noise factor)
    """
    return np.random.randn(n)*amp

def attenuate_traces(trace, fref, tstar):
    """
    Applies t* operator to obspy trace object

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
    signal = trace.data.copy()
    delta = trace.stats.delta
    nsamps = trace.stats.npts
    attenuated_signal = apply_tstar(signal, fref, tstar, delta)
    
    trace.data = attenuated_signal
    return trace
    
def apply_tstar(signal, fref, tstar, delta):
    """
    Applies a causal t* operaor (at a reference fruency fref) to the input trace

    Reference: Matheney, M. P. & Nowack, R. L. Geophys J Int 123, 1–15 (1995).
    
    Ported from msac_apply_tstar_operator by J Wookey (2022)
    
    Parameters
    ----------
    signal : TYPE
        DESCRIPTION.
    fref : TYPE
        DESCRIPTION.
    tstar : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # Take fft of trace. Supplying a larger n 0-pads trace.
    # Unlike the MATLAB implementation this is based on SciPy is 
    # not optimal for array length 2^n, instead use the Scipy next_fast_len
    # to work out the omptimal n and zero pad it (following Scipy example)
    nsamps = len(signal)
    n = next_fast_len(nsamps)
    fd_signal = fft(signal, n)
    frequencies = fftfreq(n, d=delta)[0:n//2]
    #Negative frequencies are in second half of array     
    
    #get angular frequencies
    ang_freqs = 2*np.pi*frequencies
    ang_fref = 2*np.pi*fref
    # Create causual t* multiplier 
    aw_imag = (-1j)*(1/np.pi) * tstar * ang_freqs[1:] * np.log((ang_freqs[1:])/ang_fref)
    aw_real = (-1/2)*ang_freqs[1:]*tstar
    aw = np.ones(n//2, dtype=np.complex128())
    aw[1:] = np.exp(aw_real + aw_imag)
    # Apply t* operate to frequency domain signal
    attenuated_fd_signal = np.zeros((int(n),), dtype=np.complex128())
    attenuated_pos_f = fd_signal[0:n//2]*aw
    attenuated_fd_signal[0:n//2] = attenuated_pos_f
    # Restore symettry for negative frequencies
    # index pos frequencies from 1 as we don't need to repeate f=0
    attenuated_neg_f = np.conjugate(attenuated_pos_f[1:])
    # index from [n//2+1] following scipy docs. if n is even then +/- nyquist frewquencies
    # are aliased together
    attenuated_fd_signal[n//2+1:] = np.flip(attenuated_neg_f)
    # Take inverse fft
    attenuated_signal = ifft(attenuated_fd_signal, nsamps)
    
    return attenuated_signal[0:nsamps].real

def nextpow2(i):
    n = 2
    while n < i: n = n*2
    return n

if __name__ == '__main__':
    st = obspy.read('/Users/ja17375/Projects/Matisse_Synthetics/ppv1/ideal/Noise20/data/SWAV01.BHE')
    trace = st[0]
    trace.plot()
    apply_tstar(trace, 1, 1)
    trace.plot()