#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 12:00:03 2022

Tests for synthetic.py

@author: ja17375
"""

import numpy as np
import pytest
import unittest
from numpy.testing import assert_almost_equal
from ani_attenuation.synthetics import ricker_wavelet, time_base

@pytest.mark.parametrize("delta, nsamps", [(0.1, 100), (0.5, 200), (0.2, 2000)])
def test_time_base(delta, nsamps):
    time = time_base(delta, nsamps)
    
    assert_almost_equal(np.diff(time).mean(), delta)
    assert len(time) == nsamps
    if nsamps % 2 == 1:
        assert_almost_equal(time.mean(), 0)
    elif nsamps % 2 == 0:
        time = time_base(delta, nsamps+1)
        assert_almost_equal(time.mean(), 0)
    
@pytest.mark.parametrize("nsamps, sample_rate, dfreq", [(100, 0.02, 0.1)])
def test_ricker_wavelet(nsamps, sample_rate, dfreq):
    
    test_time = sample_rate*(np.arange(0,nsamps) - nsamps//2)
    ricker = ricker_wavelet(test_time, sigma=dfreq)
    
    assert ricker[nsamps//2] == ricker.max()
    
