from labscript import *

class NIBoard(IntermediateDevice):
    # Set what types of child devices this IntermediateDevice can have:
    allowed_children = [AnalogOut, DigitalOut, AnalogIn]
    
    # Some device specific parameters:
    n_analogs = 4
    n_digitals = 32
    digital_dtype = uint32
    
    # The maximum rate that the outputs can update:
    clock_limit = 500e3
    
    # A name for the device:
    description = 'generic_NI_Board'
    
    def __init__(self, name, parent_device, clock_type, clock_terminal, MAX_name=None, acquisition_rate=0):
        # We pass the relevant parameters to the parent class's __init__ function:
        IntermediateDevice.__init__(self, name, parent_device,clock_type)
        
        # This implementation only allows analog aquisitions at a constant rate
        self.acquisition_rate = acquisition_rate
        self.clock_terminal = clock_terminal
        self.MAX_name = name if MAX_name is None else MAX_name
        self.BLACS_connection = self.MAX_name
        
    def convert_bools_to_bytes(self, digitals):
        """converts digital outputs to an array of bitfields stored
        as self.digital_dtype"""
        outputarray = [0]*self.n_digitals
        for output in digitals:
            # output.connection is the string that the user provided at
            # instantiation of the output object.  It is, by convention
            # here, port0/line<n>, where <n> is an integer from 0 to 31
            # indicating which digital output it is:
            port, line = output.connection.replace('port','').replace('line','').split('/')
            port, line  = int(port),int(line)
            if port > 0:
                raise LabscriptError('Ports > 0 on NI Boards not implemented. ' + 
                                     'Please use port 0, or file a feature request ' + 
                                     'at redmine.physics.monash.edu.au/labscript.')
            # Pack all the 1d arrays of digital output values into their appropriate spot in a list:
            outputarray[line] = output.raw_output
        # Convert this list of arrays of digital values into
        # integer bitfields (the bitfield function is located in
        # labscript.labscript)
        bits = bitfield(outputarray,dtype=self.digital_dtype)
        return bits
            
    def generate_code(self, hdf5_file):
        # By the time this function is called during compilation, most
        # of the work has already been done.  Calling the parent class's
        # generate_code method actually does nothing at the moment,
        # but this may change in the future, so you should call it anyway.
        Device.generate_code(self, hdf5_file)
        
        # Now we collect up all the output and input objects from self.child_devices:
        analogs = {}
        digitals = {}
        inputs = {}
        for device in self.child_devices:
            if isinstance(device,AnalogOut):
                analogs[device.connection] = device
            elif isinstance(device,DigitalOut):
                digitals[device.connection] = device
            elif isinstance(device,AnalogIn):
                inputs[device.connection] = device
            else:
                raise Exception('Got unexpected device.')
                
        # Now we collect up all the output.raw_output arrays from the
        # analog outputs, and load them into a numpy recarray:
        analog_out_table = empty((len(self.parent_device.times),len(analogs)), dtype=float32)
        analog_connections = analogs.keys()
        analog_connections.sort()
        analog_out_attrs = []
        for i, connection in enumerate(analog_connections):
            output = analogs[connection]
            # A bit of error checking:
            if any(output.raw_output > 10 )  or any(output.raw_output < -10 ):
                # Bounds checking:
                raise LabscriptError('%s %s '%(output.description, output.name) +
                                  'can only have values between -10 and 10 Volts, ' + 
                                  'the limit imposed by %s.'%self.name)
            # Put the 1D array of voltages into the table: 
            analog_out_table[:,i] = output.raw_output
            # Record the output terminal name to an attribute, so that
            # BLACS knows which ones to program:
            analog_out_attrs.append(self.MAX_name +'/'+connection)
            
        # Now we make a numpy recarray of all the analog input requests:
        input_connections = inputs.keys()
        input_connections.sort()
        input_attrs = []
        acquisitions = []
        for connection in input_connections:
            input_attrs.append(self.MAX_name+'/'+connection)
            for acq in inputs[connection].acquisitions:
                # Each acquisition request is a dictionary with the
                # following data, we're just putting them all in a list
                # along with the input channel they correspond to:
                acquisitions.append((connection,acq['label'],acq['start_time'],acq['end_time'],
                                     acq['wait_label'],acq['scale_factor'],acq['units']))
        # The 'a256' dtype below limits the string fields to 256
        # characters. Can't imagine this would be an issue, but to not
        # specify the string length (using dtype=str) causes the strings
        # to all come out empty.
        acquisitions_table_dtypes = [('connection','a256'), ('label','a256'), ('start',float),
                                     ('stop',float), ('wait label','a256'),('scale factor',float), ('units','a256')]
        acquisition_table= empty(len(acquisitions), dtype=acquisitions_table_dtypes)
        # OK, now we're putting them all into the numpy array:
        for i, acq in enumerate(acquisitions):
            acquisition_table[i] = acq
            
        # And finally for digital output:
        digital_out_table = []
        if digitals:
            # We convert the arrays of boolean values to a single
            # array of bitfield integers. This is how many devices need
            # their digital values programmed, though as it happens,
            # the National Instruments cards we use do not. So actually
            # this is just for storage in the HDF5 file and this process
            # is reversed when BLACS reads the data later.
            digital_out_table = self.convert_bools_to_bytes(digitals.values())
            
        # Create the required group for this device in the HDF5 file:
        grp = hdf5_file.create_group('/devices/'+self.name)
        
        # Save the analog output table, if it exists (subclasses may have zero outputs and hence an empty table):
        if all(analog_out_table.shape): # Both dimensions must be nonzero
            analog_dataset = grp.create_dataset('ANALOG_OUTS',compression=config.compression,data=analog_out_table)
            # Save the corresponding list of channels:
            grp.attrs['analog_out_channels'] = ', '.join(analog_out_attrs)
        # Save the digital output table, if it exists:
        if len(digital_out_table): # Table must be non empty
            digital_dataset = grp.create_dataset('DIGITAL_OUTS',compression=config.compression,data=digital_out_table)
            # Save the corresponding list of channels:
            grp.attrs['digital_lines'] = '/'.join((self.MAX_name,'port0','line0:%d'%(self.n_digitals-1)))
        # Save the table of acquisitions, if it exists:
        if len(acquisition_table): # Table must be non empty
            input_dataset = grp.create_dataset('ACQUISITIONS',compression=config.compression,data=acquisition_table)
            # Save the channels for analog input:
            grp.attrs['analog_in_channels'] = ', '.join(input_attrs)
            # Save the acquisition rate for analog input:
            grp.attrs['acquisition_rate'] = self.acquisition_rate
        # Save the setting for which terminal this card should expect
        # a clock input on, provided by its parent pseudoclock. BLACS
        # needs this in order to configure the device to respond to the
        # clock ticks:
        grp.attrs['clock_terminal'] = self.clock_terminal
        
