from spinapi import *
import sys
import h5py

def program_a_dds(number, group, initialamp, initialfreq):
    pb_select_dds(number)
    amps = group['AMP_REGS']
    freqs = group['FREQ_REGS']
    phases = group['PHASE_REGS']
    
    amps[0] = initialamp
    freqs[0] = initialfreq
    
    program_amp_regs(*amps)
    program_freq_regs(*freqs)
    program_phase_regs(*phases)
    
    return amps, freqs
 
 
def program_a_pulseblaster(number, group, initamp0, initfreq0, initamp1, initfreq1):
    pb_select_board(number)
    amps0, freqs0 = program_a_dds(0, group['DDS0'],initamp0, initfreq0)
    amps1, freqs1 = program_a_dds(1, group['DDS1'],initamp1, initfreq1)
    pb_start_programming(PULSE_PROGRAM)
    for args in group['PULSE_PROGRAM']:
        pb_inst_dds2(*args)
    pb_stop_programming()
    pb_close()
    #Let's get the final state of the DDS's. z's are the args I
    #don't care about:
    freqreg0,z,ampreg0,z,z,freqreg1,z,ampreg1,z,z,z,z,z,z = args
    finalfreq0 = freqs0[freqreg0]
    finalfreq1 = freqs1[freqreg1]
    finalamp0 = amps0[ampreg0]
    finalamp1 = amps1[ampreg1]
    
    return finalamp0, finalfreq0, finalamp1, finalfreq1 


if __name__== '__main__':        
    
    if not len(sys.argv) > 1 or not sys.argv[-1].endswith('.h5'):
        sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
        sys.exit(1)

    initial_states = {}
    final_states = {}
    raw_initialstates = [float(arg) for arg in sys.argv[1:-1]]
    assert len(raw_initialstates) % 4 == 0
    for i in range(len(raw_initialstates)/4):
        initial_states[i] = raw_initialstates[4*i:4*(i + 1)]
    print initial_states
    assert False
    with h5py.File(sys.argv[-1],'r') as hdf5_file:
        try:
            pb_init()
            pb_core_clock(75)
            for group in hdf5_file:
                if group.startswith('pulseblaster'):
                    number = group.split('_')[-1]
                    final_states[number] = program_a_pulseblaster(number, group, *initialstates[number])
            pb_close()
        except:
            pb_close()
            raise
    for pb in final_states:
        for item in final_states.items():
            print item,
        








