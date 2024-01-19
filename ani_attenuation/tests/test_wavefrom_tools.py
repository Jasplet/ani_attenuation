import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import unit_impulse
from ani_attenuation.waveform_tools import apply_tstar

def test_apply_tstar():

    delta_signal = unit_impulse((500,), idx='mid')
    samp_rate = 0.05 # s
    fref = 10 # Nyquist for samp rate of 0.05 s
    time_base = np.arange(0,500)*samp_rate

    delta_ts1 = apply_tstar(delta_signal, 1, samp_rate, fref)
    delta_ts2 = apply_tstar(delta_signal, 2, samp_rate, fref)
    fig = plt.figure()
    ax = fig.add_subplot(1)
    ax.plot(time_base, delta_signal)
    ax.plot(time_base, delta_ts1)
    ax.plot(time_base, delta_ts2)
    plt.show()

if __name__ == '__main__':
    test_apply_tstar()