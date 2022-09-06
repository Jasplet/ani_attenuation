#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 09:51:11 2022

@author: ja17375
"""

import numpy as np
import obspy

from scipy.interpolate import PchipInterpolator

def rotate_traces(trace1, trace2, theta):
    '''
    Rotate two obspy traces to the desired rotation angle theta. 
    
    Adapted port of msac_rotate 
    Parameters
    ----------
    trace1 : obspy Trace
        y-direction trace (e.g. North component)
    trace2 : obspy Trace
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
    ydata = trace1.data.copy()
    xdata = trace2.data.copy()
    trace1p = trace1.copy()
    trace2p = trace2.copy()
    
    ydatap, xdatap = rotate(ydata, xdata, theta)
    trace1p.data = ydatap
    trace2p.data = xdatap
    trace1p.stats.sac['cmpaz'] = trace1.stats.sac['cmpaz'] + theta 
    trace2p.stats.sac['cmpaz'] = trace2.stats.sac['cmpaz'] + theta
    return trace1p, trace2p
    
def rotate(x, y, theta):
    """
    Rotate 2 orthogonal traces

    Parameters
    ----------
    tr1 : TYPE
        DESCRIPTION.
    tr2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if len(y) != len(x):
        raise ValueError('Traces must be same length!')
    if (len(y.shape) > 1 ) or (len(x.shape) > 1):
        raise ValueError('Traces must be 1-D arrays!')
    
    trs = np.vstack((x, y)) # tr1 is y so goes second
    # form rotation matrix
    radtheta = np.deg2rad(theta)
    c = np.cos(radtheta)
    s = np.sin(radtheta)
    rotmat = np.array([
        [ np.cos(radtheta), np.sin(radtheta)],
        [-np.sin(radtheta), np.cos(radtheta)]
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
    attenuated_signal = apply_tstar(signal, fref, tstar, delta, nsamps)
    
    trace.data = attenuated_signal
    return trace
    
def apply_tstar(signal, fref, tstar, delta, nsamps):
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
    n = int(nextpow2(nsamps))
    fd_signal = np.fft.fft(signal, n)
    frequencies = np.fft.fftfreq(n, d=delta)[0:n//2]
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
    
    return attenuated_signal[0:nsamps].real

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