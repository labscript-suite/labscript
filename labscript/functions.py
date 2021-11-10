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
    """Print time with a descriptive string.

    Useful debug tool to print time at a specific point
    in the shot, during shot compilation. Helpful when
    the time is calculated.

    Args:
        t (float): Time to print
        description (str): Descriptive label to print with it
    """
    print('t = {0:.9f} s:'.format(t),description)

def ramp(duration, initial, final):
    """Defines a linear ramp.

    f(t) = (final - initial)*t/duration + initial

    Args:
        duration (float): Duration of ramp
        initial (float): Starting value of ramp
        final (float): Ending value of ramp

    Returns:
        func: Function that takes a single parameter `t`.
    """
    m = (final - initial)/duration
    return lambda t: m*t + initial

def sine(duration, amplitude, angfreq, phase, dc_offset):
    """Defines a sine wave.

    f(t) = amplitude*sin(angfreq*t + phase) + dc_offset

    Args:
        duration (float): Not used.
        amplitude (float): Amplitude of sine wave.
        angfreq (float): Angular frequency of sine wave.
        phase (float): Phase of sine wave.
        dc_offset (float): Verticle offset of sine wave.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    return lambda t: amplitude*sin(angfreq*(t) + phase) + dc_offset
    
def sine_ramp(duration, initial, final):
    """Defines a square sinusoidally increasing ramp.

    f(t) = (final-initial)*(sin(pi*t/(2*duration)))^2 + initial

    Args:
        duration (float): Length of time for the ramp to complete.
        initial (float): Initial value of ramp.
        final (float): Final value of ramp.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    return lambda t: (final-initial)*(sin(pi*(t)/(2*duration)))**2 + initial
    
def sine4_ramp(duration, initial, final):
    """Defines a quartic sinusoidally increasing ramp.

    f(t) = (final-initial)*(sin(pi*t/(2*duration)))^4 + initial

    Args:
        duration (float): Length of time for the ramp to complete.
        initial (float): Initial value of ramp.
        final (float): Final value of ramp.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    return lambda t: (final-initial)*(sin(pi*(t)/(2*duration)))**4 + initial
    
def sine4_reverse_ramp(duration, initial, final):
    """Defines a quartic sinusoidally decreasing ramp.

    f(t) = (final-initial)*(sin(pi/2+pi*t/(2*duration)))^4 + initial

    Args:
        duration (float): Length of time for the ramp to complete.
        initial (float): Initial value of ramp.
        final (float): Final value of ramp.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    return lambda t: (final-initial)*(sin(pi/2+pi*(t)/(2*duration)))**4 + initial
    
def exp_ramp(duration,initial,final,zero):
    """Defines an exponential ramp via offset value.

    f(t) = (initial-zero)*e^(-rate*t) + zero
    rate = log((initial-zero)/(final-zero))/duration

    Args:
        duration (float): Length of time for the ramp to complete
        initial (float): Initial value of ramp.
        final (float): Final value of ramp.
        zero (float): Zero offset of ramp.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    rate = 1/duration * log((initial-zero)/(final-zero))
    return lambda t: (initial-zero)*exp(-rate*(t)) + zero
    
def exp_ramp_t(duration,initial,final,time_constant):
    """Defines an exponential ramp via time constant.

    f(t) = (initial-zero)*e^(-t/time_constant) + zero
    zero = (final-initial*e^(-duration/time_constant))/(1-e^(-duration/time_constant))

    Args:
        duration (float): Length of time for the ramp to complete
        initial (float): Initial value of ramp.
        final (float): Final value of ramp.
        zero (float): Zero offset of ramp.

    Returns:
        func: Function that takes a single parameter `t`.
    """
    zero = (final-initial*exp(-duration/time_constant)) / (1-exp(-duration/time_constant))
    return lambda t: (initial-zero)*exp(-(t)/time_constant) + zero

def piecewise_accel(duration,initial,final):
    """Defines a piecewise acceleration.

    Args:
        duration (float): Length of time for the acceleration to complete.
        initial (float): Initial value.
        final (float): Final value.
    """
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
    """Returns a function that interpolates a pulse sequence.

    Relies on :obj:`numpy.digitize` to perform the interpolation.

    Args:
        pulse_sequence (:obj:`numpy:numpy.ndarray`): 2-D timeseries of
            change times and associated states.
        period (float): How long, in seconds, to hold the final state
            before repeating the sequence.

    Returns:
        func: Interpolating function that takes a single parameter `t`.
        Only well defined if `t` falls within the `pulse_sequence` change times.
    """
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
