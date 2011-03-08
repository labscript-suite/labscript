from spinapi import *
import sys
import h5py

if not len(sys.argv) > 1:
    sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
    sys.exit(1)
    
with h5py.File(sys.argv[-1],'r') as hdf5_file:
    try:
        pb_init()
        pb_core_clock(75)
        pb_select_dds(0)
        pb_start_programming(PULSE_PROGRAM)
        for args in hdf5_file['/PulseBlaster/PULSE_PROGRAM']:
            pb_inst_dds2(*args)
        pb_stop_programming()
        pb_close()
    except:
        pb_close()
        raise
