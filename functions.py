#####################################################################
#                                                                   #
# /functions.py                                                     #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

from __future__ import division
from pylab import *

def ramp(duration, initial, final):
    m = (final - initial)/duration
    return lambda t: m*t + initial

def sine(duration, amplitude, angfreq, phase, dc_offset):
    return lambda t: amplitude*sin(angfreq*(t) + phase) + dc_offset
    
def sine_ramp(duration, initial, final):
    return lambda t: (final-initial)*(sin(pi*(t)/(2*duration)))**2 + initial
    
def sine4_ramp(duration, initial, final):
    return lambda t: (final-initial)*(sin(pi*(t)/(2*duration)))**4 + initial
    
def sine4_reverse_ramp(duration, initial, final):
    return lambda t: (final-initial)*(sin(pi/2+pi*(t)/(2*duration)))**4 + initial
    
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
    
def pulse_sequence(pulse_sequence,period):    
    def pulse_function(t):
        try:
            len(t)
            is_array = True
        except TypeError:
            t = array([t])
            is_array = False
            
        times = t%period
        states = zeros(len(t),dtype = int)
        for i, a_time in enumerate(times):
            for time, state in sorted(pulse_sequence, key=lambda x: x[0], reverse=True):
                if a_time >= time:
                    states[i] = state
                    break
                    
        if is_array:
            return states
        else:
            return states[0]
    return pulse_function
