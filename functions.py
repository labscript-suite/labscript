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
    
def exp_ramp(t,duration,initial,final,zero):
    rate = 1/duration * log((initial-zero)/(final-zero))
    return lambda x: (initial-zero)*exp(-rate*(x-t)) + zero
    
def exp_ramp_t(t,duration,initial,final,time_constant):
    zero = (final-initial*exp(-duration/time_constant)) / (1-exp(-duration/time_constant))
    return lambda x: (initial-zero)*exp(-(x-t)/time_constant) + zero
