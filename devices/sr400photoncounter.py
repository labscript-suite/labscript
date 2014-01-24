#####################################################################
#                                                                   #
# /devices/sr400photoncounter.py                                    #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

from labscript import DigitalOut, LabscriptError

class SR400PhotonCounter(DigitalOut):
    description = 'SRS400 Gated Photon Counter'
    def __init__(self, name, parent_device, connection, com_port):
        DigitalOut.__init__(self, name, parent_device, connection)
        self.BLACS_connection = com_port
        self.duration = None
        
    def acquire(self, t, bin_size, n_periods, dwell_time=2e-3):
        if self.duration is not None:
            raise LabscriptError('already triggered')
        self.t = t
        self.bin_size_10MHz = bin_size*1e7
        self.n_periods = int(n_periods)
        assert dwell_time >= 2e-3
        self.dwell_time = dwell_time
        self.go_high(t)
        # A short pulse. Potential problems, feel free to change length:
        self.go_low(t+1e-3)
        duration = (bin_size + dwell_time)*n_periods
        self.duration = duration
        return duration
        
    def generate_code(self, hdf5_file):
        if self.duration is None:
            return
        if self.t + self.duration > self.pseudoclock.stop_time:
            raise LabscriptError('%s acquires past the end of the experiment!'%self.name)
            
        group = hdf5_file.create_group('/devices/%s'%self.name)
        group.attrs['bin_size_10MHz'] = self.bin_size_10MHz
        group.attrs['n_periods'] = self.n_periods
        group.attrs['dwell_time'] = self.dwell_time

