from __future__ import division
from pylab import *

def ramp(t, duration, initial, final):
    m = (final - initial)/duration
    c = initial - m*t
    return lambda x: m*x + c

def sine(t, duration, amplitude, angfreq, phase, dc_offset):
    return lambda x: amplitude*sin(angfreq*(x-t) + phase) + dc_offset
    
def sine_ramp(t, duration, initial, final):
    return lambda x: (final-initial)*(sin(pi*(x-t)/(2*duration)))**2 + initial
