from __future__ import division
from pylab import *

def ramp(duration, initial, final):
    m = (final - initial)/duration
    return lambda t: m*t + initial

def sine(duration, amplitude, angfreq, phase, dc_offset):
    return lambda t: amplitude*sin(angfreq*(t) + phase) + dc_offset
    
def sine_ramp(duration, initial, final):
    return lambda t: (final-initial)*(sin(pi*(t)/(2*duration)))**2 + initial
    
def exp_ramp(duration,initial,final,zero):
    rate = 1/duration * log((initial-zero)/(final-zero))
    return lambda t: (initial-zero)*exp(-rate*(t)) + zero
    
def exp_ramp_t(duration,initial,final,time_constant):
    zero = (final-initial*exp(-duration/time_constant)) / (1-exp(-duration/time_constant))
    return lambda t: (initial-zero)*exp(-(t)/time_constant) + zero


def piecewise_accel(duration,initial,final):
    a = (final-initial)
    return lambda t: initial + a * (
    (9./2 * t**3/duration**3) * (t<duration/3)
    + (-9*t**3/duration**3 + 27./2*t**2/duration**2 - 9./2*t/duration + 1./2) * (t<2*duration/3)*(t>=duration/3)
    + (9./2*t**3/duration**3 - 27./2 * t**2/duration**2 + 27./2*t/duration - 7./2) * (t>= 2*duration/3))