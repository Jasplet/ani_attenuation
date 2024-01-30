#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 10:32:28 2022

@author: ja17375
"""

import numpy as np
import obspy

from ani_attenuation.waveform_tools import rotate, apply_tshift, time_base, randn_noise, attenuate_traces, rotate_traces
    
def gen_synthetic_split(fast, tlag, **kwargs):
    '''
    Function to generate a synthetic shear-wave

    Returns
    -------
    None.
    '''
    defaults = {'noise':0, 'dfreq':0.2, 'delta':0.05, 'spol':30}
    if 'noise'in kwargs:
        noise = kwargs['noise']
    else:
        noise = defaults['noise']
    
    if 'dfreq' in kwargs:
        dfreq = kwargs['dfreq']
    else:
        dfreq = defaults['dfreq']
    
    if 'delta' in kwargs:
        delta = kwargs['delta']
    else:
        delta = defaults['delta']
    if 'spol' in kwargs:
        spol = kwargs['spol']
    else:
        spol = defaults['spol']
    if 'fref' in kwargs:
        fref = kwargs['fref']
    else:
        # Set fref to Nyquist frequency
        fref = 1/(2*delta)
    if 'nsamps' in kwargs:
        nsamps = kwargs['nsamps']
    else:
        #default number of samples is 10 times dominant preiod (1/dfreq)
        nsamps = int(1 + 10*(1/ dfreq) / delta)

    time = time_base(delta, nsamps)
    # create wavelet
    wavelet = gabor_wavelet(time, dfreq)
    max_amp = np.max(np.abs(wavelet))
    waveletN = wavelet*np.cos(np.deg2rad(spol)) + randn_noise(int(nsamps), max_amp*noise)
    waveletE = wavelet*np.sin(np.deg2rad(spol)) + randn_noise(int(nsamps), max_amp*noise)
    waveletZ = randn_noise(int(nsamps), max_amp*noise)
    # Take N, E Wavelets and add splitting    
    waveletF, waveletS = rotate(waveletN, waveletE, fast)
    waveletS = apply_tshift(waveletS, tlag, delta, time[0])
    waveletN, waveletE = rotate(waveletF, waveletS, -1*fast)
   
    # Now add metadata needed to make a obspy Trace/Stream object
    stats = make_stats_dict(delta, nsamps, dfreq, time, spol)
    traceN = obspy.Trace()
    traceN.stats = stats.copy()
    traceN.stats.channel = 'BHN'
    traceN.stats.sac.cmpaz = 0
    traceN.stats.sac.cmpinc = 90
    traceN.stats.sac.kcmpnm = 'BHN'
    traceN.data = waveletN
    #East cmp
    traceE = obspy.Trace()
    traceE.stats = stats.copy()
    traceE.stats.channel = 'BHE'
    traceE.stats.sac.cmpaz = 90
    traceE.stats.sac.cmpinc = 90
    traceE.stats.sac.kcmpnm = 'BHE'
    traceE.data = waveletE
    # Vertical cmp
    traceZ = obspy.Trace()
    traceZ.stats = stats.copy()
    traceZ.stats.channel = 'BHZ'
    traceZ.stats.sac.cmpaz = 00
    traceZ.stats.sac.cmpinc = 0
    traceZ.stats.sac.kcmpnm = 'BHZ'
    traceZ.data = waveletZ
    # attenuate traces if needed 
    if 'dtstar' in kwargs:
        dtstar = kwargs['dtstar']
        [traceF, traceS] = rotate_traces(traceN, traceE, fast)
        if dtstar > 0:
            traceSA = attenuate_traces(traceS, dtstar, fref)
            [traceN, traceE] = rotate_traces(traceF, traceSA, -1*fast)
        elif kwargs['dtstar'] < 0:
            traceFA = attenuate_traces(traceF, abs(kwargs['dtstar']),fref)
            [traceN, traceE] = rotate_traces(traceFA, traceS, -1*fast)
        elif dtstar == 0:
            print('dt* = 0')
            [traceN, traceE] = rotate_traces(traceF, traceS, -1*fast)
        else:
            raise ValueError(f'Unknown dt* {dtstar}')
    synthetic = obspy.Stream([traceN, traceE, traceZ])
    return synthetic

def gabor_wavelet(t, dfreq, gamma=6, v=np.pi*(2/5), t0=0):
    '''
    Port of gabor wavelet function from msac_splitwave2
    '''
    term1 = 2*np.cos(2. * np.pi * dfreq * (t - t0) + v)
    term2 = np.exp((-4 * (np.pi**2) * (dfreq**2) * (t-t0)**2)/gamma)
    wavelet = term1*term2
    
    return wavelet / max(wavelet)
    
def ricker_wavelet(t, sigma=0.1):
    '''
    Parameters
    ----------
    t : array
        Time base for wavelet 
    sigma :float, optional
        dominant frequency of wavelet. The default is 0.1 [Hz]

    Returns
    -------
    ricker : array
        ricker wavelet

    '''
    term1 = 2 / (np.sqrt(3*sigma)*(np.pi**0.25))
    term2 = 1 - (t / sigma)**2
    term3 = np.exp(-0.5 * (t / sigma)**2)
    ricker =  term1*term2*term3
    
    return ricker

def make_stats_dict(delta, nsamps, dfreq, time, spol):
    '''
    Makes base stats object for obspy trace
    '''
    stats = obspy.core.trace.Stats()
    sachdrs = obspy.core.util.attribdict.AttribDict()
    stats.station = 'SYN'
    stats.starttime = obspy.UTCDateTime('3001-01-01T00:00:00')
    stats._format = 'SAC'
    stats.sampling_rate = 1/delta
    stats.delta = delta
    sachdrs.delta = delta
    sachdrs.b = time[0]
    sachdrs.e = time[-1]
    sachdrs.baz = spol # Imposing assumption that BAZ == SPOL as in SKS case
    sachdrs.evla = 0
    sachdrs.evlo = 0
    sachdrs.evdp = 500
    sachdrs.stlo = 0
    sachdrs.stla = 80
    sachdrs.kstnm = 'SYN'
    #Set defualt windows (double dominant period 1/dfreq)
    sachdrs.a = -1*(1/dfreq)
    sachdrs.f = (1/dfreq)
    sachdrs.user0 = sachdrs.a
    sachdrs.user2 = sachdrs.f
    stats.sac = sachdrs
    return stats

if __name__ == '__main__':
    wv = gen_synthetic_split(45, 1, spol=20, dtstar=1.1, dfreq=0.1)
    wv.plot()