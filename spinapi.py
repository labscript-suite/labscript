import ctypes
import time
from spinconfig import config

spinapi = ctypes.cdll.LoadLibrary(config['spinapi'])

# Defines for different pb_inst instruction types
CONTINUE = 0
STOP = 1
LOOP = 2
END_LOOP = 3
JSR = 4
RTS = 5
BRANCH = 6
LONG_DELAY = 7
WAIT = 8
RTI = 9

# Defines for using different units of time
ns = 1.0
us = 1000.0
ms = 1000000.0
s  = 1000000000.0

# Defines for using different units of frequency
MHz = 1.0
kHz = .001
Hz = .000001

# Defines for start_programming
PULSE_PROGRAM  = 0
FREQ_REGS = 1
PHASE_REGS = 2

# Defines for enabling analog output
ANALOG_ON = 1
ANALOG_OFF = 0

# Defines for resetting the phase:
PHASE_RESET = 1
NO_PHASE_RESET = 0

def spinpts_get_version():
    spinapi.spinpts_get_version.restype = ctypes.c_char_p
    return spinapi.spinpts_get_version()
    
def pb_get_error():
    spinapi.pb_get_error.restype = ctypes.c_char_p
    return spinapi.pb_get_error()

def pb_status_message():
    spinapi.pb_status_message.restype = ctypes.c_char_p
    message = spinapi.pb_status_message()
    return message
    
def pb_count_boards():
    spinapi.pb_count_boards.restype = ctypes.c_int
    result = spinapi.pb_count_boards()
    if result == -1: raise RuntimeError(pb_get_error())  
    return result

def pb_select_board(board_num):
    spinapi.pb_select_board.restype = ctypes.c_int
    retcode = spinapi.pb_select_board(ctypes.c_int(board_num))
    if retcode != 0: raise RuntimeError(pb_get_error())
                               
def pb_init():
    spinapi.pb_init.restype = ctypes.c_int
    retcode = spinapi.pb_init()
    if retcode != 0: raise RuntimeError(pb_get_error())

def pb_core_clock(clock_freq):
    spinapi.pb_core_clock.restype = ctypes.c_int
    retcode = spinapi.pb_core_clock(ctypes.c_double(clock_freq))
    if retcode != 0: raise RuntimeError(pb_get_error())
                       
def pb_start_programming(device):
    spinapi.pb_start_programming.restype = ctypes.c_int
    retcode = spinapi.pb_start_programming(ctypes.c_int(device))
    if retcode != 0: raise RuntimeError(pb_get_error())

def pb_select_dds(dds):
    spinapi.pb_select_dds.restype = ctypes.c_int
    result = spinapi.pb_select_dds(ctypes.c_int(dds))
    if result < 0: raise RuntimeError(pb_get_error())
    
def pb_set_phase(phase):
    spinapi.pb_set_phase.restype = ctypes.c_int
    result = spinapi.pb_set_phase(ctypes.c_double(phase))
    if result < 0: raise RuntimeError(pb_get_error())

def pb_set_freq(freq):
    spinapi.pb_set_freq.restype = ctypes.c_int
    result = spinapi.pb_set_freq(ctypes.c_double(freq))
    if result < 0: raise RuntimeError(pb_get_error())

def pb_set_amp(amp, register):
    spinapi.pb_set_amp.restype = ctypes.c_int
    result = spinapi.pb_set_amp(ctypes.c_float(amp),ctypes.c_int(register))
    if result < 0: raise RuntimeError(pb_get_error())
    
def pb_inst(flags, inst, inst_data, length):
    """This is a convenience function, with all the analogue stuff set to zero.
       Some of the documentation refers to a function with this name and syntax,
       even though it's not actually in the API. It's modified here to take a 
       string of ones and zeros for the 'flags' argument, converting it to an 
       int for the C function call."""
    spinapi.pb_inst_dds2.restype = ctypes.c_int
    z = ctypes.c_int(0)
    # [::-1] reverses any iterable. This is Python idiom...despite not being obvious at all.
    result = spinapi.pb_inst_dds2(z,z,z,z,z,z,z,z,z,z, ctypes.c_int(int(flags[::-1],2)),
                                  ctypes.c_int(inst),ctypes.c_int(inst_data),ctypes.c_double(length))
    if result < 0: raise RuntimeError(pb_get_error())
    return result

def pb_inst_dds2(freq0,phase0,amp0,dds_en0,phase_reset0,
                 freq1,phase1,amp1,dds_en1,phase_reset1,
                 flags, inst, inst_data, length):
    """Gives a full instruction to the pulseblaster, with DDS included. The flags argument can be
       either an int representing the bitfield for the flag states, or a string of ones and zeros."""
    spinapi.pb_inst_dds2.restype = ctypes.c_int
    if isinstance(flags, str):
        flags = int(flags[::-1],2)
    result = spinapi.pb_inst_dds2(ctypes.c_int(freq0),ctypes.c_int(phase0),ctypes.c_int(amp0),
                                  ctypes.c_int(dds_en0),ctypes.c_int(phase_reset0),
                                  ctypes.c_int(freq1),ctypes.c_int(phase1),ctypes.c_int(amp1),
                                  ctypes.c_int(dds_en1),ctypes.c_int(phase_reset1),
                                  ctypes.c_int(),ctypes.c_int(inst),
                                  ctypes.c_int(inst_data),ctypes.c_double(length))
    if result < 0: raise RuntimeError(pb_get_error())
    return result

# More convenience functions:
def program_freq_regs(*freqs):
    pb_start_programming(FREQ_REGS)
    for freq in freqs:
        pb_set_freq(freq)
    pb_stop_programming()
    if len(freqs) == 1:
        return 0
    else:
        return tuple(range(len(freqs)))

def program_phase_regs(*phases):
    pb_start_programming(PHASE_REGS)
    for phase in phases:
        pb_set_phase(phase)
    pb_stop_programming()
    if len(phases) == 1:
        return 0
    else:
        return tuple(range(len(phases)))

def program_amp_regs(*amps):
    for i, amp in enumerate(amps):
        pb_set_amp(amp,i)
    if len(amps) == 1:
        return 0
    else:
        return tuple(range(len(amps)))

def pb_stop_programming():
    spinapi.pb_stop_programming.restype = ctypes.c_int
    retcode = spinapi.pb_start_programming()
    if retcode != 0: raise RuntimeError(pb_get_error())

def pb_start():
    spinapi.pb_start.restype = ctypes.c_int
    retcode = spinapi.pb_start()
    if retcode != 0: raise RuntimeError(pb_get_error())

def pb_stop():
    spinapi.pb_stop.restype = ctypes.c_int
    retcode = spinapi.pb_stop()
    if retcode != 0: raise RuntimeError(pb_get_error())
                                                        
def pb_close():
    spinapi.pb_close.restype = ctypes.c_int
    retcode = spinapi.pb_close()
    if retcode != 0: raise RuntimeError(pb_get_error())

def program_and_run(program):
    """This bit has nothing to do with the C API and is for programming the PulseBlaster
    with a Python script that imports this one. See the if __name__ == '__main__' block below
    for an example."""
    print 'Spincore API version', spinpts_get_version()
    print 'Initialising board...',
    pb_init()
    pb_core_clock(50)
    print 'done.'
    print 'Programming board...',
    program()
    print 'done.'
    print 'Starting program...',
    pb_start()
    print 'press Ctrl-C to stop.'
    print pb_status_message().replace('\n','')
    try:
        while pb_status_message() == 'Board is running.\n':
            time.sleep(0.5)
        print pb_status_message().replace('\n','')
    except KeyboardInterrupt:
        print 'Stopping board ...',
    pb_stop()
    pb_close()
    print 'done.'

if __name__ == '__main__':

    # Use 'from spinapi import *' to import this module to another script in the same directory.
    # Then you can define your program like so:
    
    def program():
        pb_start_programming(PULSE_PROGRAM)
        start = pb_inst(    '111111111111',    CONTINUE,      0,    100.0*us    )
        _     = pb_inst(    '000000000000',    BRANCH,    start,    100.0*us    )
        pb_stop_programming()
    # And run it thusly:
    program_and_run(program)










    
