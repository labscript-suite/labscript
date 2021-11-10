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

from pylab import *
import numpy as np

def print_time(t, description):
    print('t = {0:.9f} s:'.format(t),description)

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

def square_wave(duration, level_0, level_1, frequency, phase, duty_cycle):
    def square_wave_fixed_parameters(t):
        # Phase goes from 0 to 1 (NOT 2 pi) over one period.
        edge_phase_0_to_1 = duty_cycle
        wrapped_phases = (frequency * t + phase) % 1.0
        # Ensure wrapped_phases is an array.
        wrapped_phases = np.array(wrapped_phases)

        # Round phases to avoid issues with numerics. Rounding the phase only
        # changes the output when the phase is just below a threshold where the
        # output changes values. So if a phase is just below the threshold where
        # the output changes state (within PHASE_TOLERANCE), round it up so that
        # the output does change state there. The value of PHASE_TOLERANCE is
        # based on the fact that labscript internally rounds all times to
        # multiples of 0.1 ns.
        LABSCRIPT_TIME_RESOLUTION = 0.1e-9  # 0.1 ns.
        MIN_PHASE_STEP = frequency * LABSCRIPT_TIME_RESOLUTION
        PHASE_TOLERANCE = MIN_PHASE_STEP / 2.0
        # Round phases near level_0 -> level_1 transition at phase =
        # edge_phase_0_to_1.
        is_near_edge = np.isclose(
            wrapped_phases,
            edge_phase_0_to_1,
            rtol=0,
            atol=PHASE_TOLERANCE,
        )
        wrapped_phases[is_near_edge] = edge_phase_0_to_1
        # Round phases near level_1 -> level_0 transition at phase = 1.
        is_near_edge = np.isclose(
            wrapped_phases,
            1,
            rtol=0,
            atol=PHASE_TOLERANCE,
        )
        wrapped_phases[is_near_edge] = 0

        # Initialize array to store output values.
        outputs = np.full_like(t, level_0)

        # Use boolean indexing to set output to level_1 at the appropriate
        # times. For example level_0 for phases [0, 0.5) and level_1 for phases
        # [0.5, 1.0) when duty_cycle is 0.5.
        level_1_times = (wrapped_phases >= edge_phase_0_to_1)
        outputs[level_1_times] = level_1
        return outputs
    return square_wave_fixed_parameters

def pulse_sequence(pulse_sequence,period):
    pulse_sequence = np.asarray(sorted(pulse_sequence, key=lambda x: x[0], reverse=True))
    pulse_sequence_times = pulse_sequence[:, 0]
    pulse_sequence_states = pulse_sequence[:, 1]

    def pulse_function(t):
        try:
            len(t)
            is_array = True
        except TypeError:
            t = array([t])
            is_array = False

        times = t % period
        indices = np.digitize(times, pulse_sequence_times, right=False)
        states = pulse_sequence_states[indices]

        if is_array:
            return states
        else:
            return states[0]

    return pulse_function
