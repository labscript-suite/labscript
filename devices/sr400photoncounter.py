from labscript import DigitalOut

class SR400PhotonCounter(DigitalOut):
    description = 'SRS400 Gated Photon Counter'
    def __init__(self, name, parent_device, connection, com_port):
        DigitalOut.__init__(self, name, parent_device, connection)
        self.BLACS_connection = com_port
        self.duration = None
        
    def acquire(self, t, duration):
        if self.duration is not None:
            raise LabscriptError('already triggered')
        self.duration = duration
        self.go_high(t)
        self.go_low(t + duration)
        return duration
        
    def generate_code(self, hdf5_file):
        if self.duration is None:
            raise LabscriptError('was never triggered')
        group = hdf5_file.create_group('/devices/%s'%self.name)
        group.attrs['duration'] = self.duration
        
