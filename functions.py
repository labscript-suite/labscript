from pylab import *

def ramp(t,duration,initial,final):
    m = (final - initial)/float(duration) # be sure to prevent integer division!
    c = initial - m*t
    return lambda x: m*x + c

def sine(t,amplitude,angfreq,phase,dc_offset):
    return lambda x: amplitude*sin(angfreq*(x-t) + phase) + dc_offset
