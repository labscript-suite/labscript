import os
import sys
import keyword

import h5py
from pylab import *

import functions


def fastflatten(inarray):
    """A faster way of flattening our arrays than pylab.flatten.
    pylab.flatten returns a generator which takes a lot of time and memory
    to convert into a numpy array via array(list(generator)).  The problem
    is that generators don't know how many values they'll return until
    they're done. This algorithm produces a numpy array directly by
    first calculating what the length will be. It is several orders of
    magnitude faster. Note that we can't use numpy.ndarray.flatten here
    since our inarray is really a list of 1D arrays of varying length
    and/or single values, not a N-dimenional block of homogeneous data
    like a numpy array."""
    total_points = sum([len(element) if iterable(element) else 1 for element in inarray])
    flat = empty(total_points,dtype=float32)
    i = 0
    for val in inarray:
        if iterable(val):
            flat[i:i+len(val)] = val[:]
            i += len(val)
        else:
            flat[i] = val
            i += 1
    return flat
    
def discretise(t,y,stop_time):
    tnew = zeros((len(t),2))
    ynew = zeros((len(y),2))

    tnew[:,0] = t[:]
    tnew[:-1,1] = t[1:]
    tnew= tnew.flatten()
    tnew[-1] = stop_time

    ynew[:,0] = y[:]
    ynew[:,1] = y[:]
    ynew= ynew.flatten()[:]
    return tnew, ynew

def plot_outputs(devices='all',display=False):
    if devices == 'all':
        devices = inventory
    for device in devices:
        for i, output in enumerate(device.outputs):
            t,y = discretise(device.flat_times,output.raw_output, device.stop_time)
            plot(t,y,'-',label=output.name)

    grid(True)
    xlabel('time (seconds)')
    ylabel('output values')
    title('Pseudoclocked outputs')
    legend(loc='upper left')
    axis([0,max([device.stop_time for device in inventory]),-1,5.5])
    if display:
        show()
    else:
        savefig('outputs.png')
    
    
class IODevice:
    """This class represents a grouping of outputs (analogue or
    digital) that all share a common pseudoclock. When
    IODevice.make_instruction_table is called, instructions are collected
    from the IODevice's attached outputs, and arrays of values are
    generated for them to output on each clock tick, as well as
    instructions for the pseudo clock itself."""

    # Maximum clock rate of the pseudo clock. To be overridden by subclasses:
    clock_limit = 1e9
    description = 'IO device'
    
    def __init__(self,name,stop_time=20):
        self.outputs = []
        self.name = name
        self.stop_time = stop_time
        inventory.append(self)
        
    def collect_change_times(self):
        """Asks all connected outputs for a list of times that they
        change state. Takes the union of all of these times. Note
        that at this point, a change from holding-a-constant-value
        to ramping-through-values is considered a single state
        change. The clocking times will be filled in later in the
        expand_change_times function, and the ramp values filled in with
        expand_timeseries."""
        # Use a set to avoid duplicates:
        self.change_times = set()
        for output in self.outputs:
            output.make_times()
            self.change_times = self.change_times.union(output.times)
        self.change_times = list(self.change_times)
        self.change_times.sort()
        # Check that no two instructions are too close together:
        for i, t in enumerate(self.change_times[:-1]):
            dt = self.change_times[i+1] - t
            if dt < 1.0/self.clock_limit:
                sys.stderr.write('ERROR: Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(self.change_times[i+1])) +
                                  'One or more connected devices cannot support update delays shorter than %s sec. Stopping.\n'%str(1.0/self.clock_limit))
                sys.exit(1)
        
    def make_timeseries(self):
        """Instructs each connected output to construct a list of its
        states at each time in self.change_times, that is at any point in
        time that one or more connected outputs change state. By state,
        I don't mean the value of the output at that moment, rather
        I mean what instruction it has. This might be a single value,
        or it might be a reference to a function for a ramp etc."""
        for output in self.outputs:
            output.make_timeseries(self.change_times)
        
    def expand_change_times(self):
        """For each time interval delimited by self.change_times,
        constructs an array of times at which the clock for this device
        needs to tick. If the interval has all outputs having constant
        values, then only the start time is stored.  If one or more
        outputs are ramping, then the clock ticks at the maximum clock
        rate requested by any of the outputs. Also produces a higher
        level description of the clocking; self.clock. This list contains
        the information that facilitates programming a pseudo clock
        using loops."""
        self.all_times = []
        self.clock = []
        for i, time in enumerate(self.change_times):
            # what's the fastest clock rate?
            maxrate = 0
            for output in self.outputs:
                # Check if output is sweeping and has highest clock rate
                # so far. If so, store its clock rate to max_rate:
                if isinstance(output.timeseries[i],dict) and output.timeseries[i]['clock rate'] > maxrate:
                    # It does have the highest clock rate? Then store that rate to max_rate:
                    maxrate = output.timeseries[i]['clock rate']
            if maxrate:
                # If there was ramping at this timestep, how many clock ticks fit before the next instruction?
                n_ticks, remainder = divmod((self.change_times[i+1] - time)*maxrate,1)
                n_ticks = int(n_ticks)
                # Can we squeeze the final clock cycle in at the end?
                if remainder and remainder/float(maxrate) >= 1/float(self.clock_limit):
                    # Yes we can. Clock speed will be as
                    # requested. Otherwise the final clock cycle will
                    # be too long, by the fraction 'remainder'.
                    n_ticks += 1
                duration = n_ticks/float(maxrate) # avoiding integer division
                ticks = linspace(time,time + duration,n_ticks,endpoint=False)
                self.all_times.append(array(ticks,dtype=float32))
                # Note that even though all arrays are single precision
                # floating point, the numbers stored to the clock list
                # below are double precision. This is important so that
                # rounding errors in the stepsize don't become significant
                # after many clock cycles.
                if n_ticks > 1:
                    # If n_ticks is only one, then this step doesn't do
                    # anything, it has reps=0. So we should only include
                    # it if n_ticks > 1:
                    self.clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate)})
                # The last clock tick has a different duration depending on the next step:
                self.clock.append({'start': ticks[-1], 'reps': 1, 'step': self.change_times[i+1] - ticks[-1]})
            else:
                self.all_times.append(time)
                try: 
                    # If there was no ramping, here is a single clock tick:
                    self.clock.append({'start': time, 'reps': 1, 'step': self.change_times[i+1] - time})
                except IndexError:
                    if i != len(self.change_times) - 1:
                        raise
                    # There is no next instruction. Hold the last clock
                    # tick until self.stop_time.
                    if self.stop_time > time:
                        self.clock.append({'start': time, 'reps': 1, 'step': self.stop_time - time})
                    # Error if self.stop_time has been set to less
                    # than the time of the last instruction:
                    elif self.stop_time < time:
                        sys.stderr.write('ERROR: %s %s has had its stop time set to earlier than its last instruction!'%(self.description,self.name))
                        sys.exit(1)
                    # If self.stop_time is the same as the time of the last
                    # instruction, then we'll get the last instruction
                    # out still, so that the total number of clock
                    # ticks matches the number of data points in the
                    # Output.raw_output arrays. We'll make this last
                    # cycle be at ten times the maximum step duration.
                    else:
                        self.clock.append({'start': time, 'reps': 1, 'step': 10.0/self.clock_limit})
                        
                        
    def expand_timeseries(self):
        for output in self.outputs:
            output.make_outputarray(self.all_times)
            
    def make_raw_output(self):
        self.flat_times = fastflatten(self.all_times)
        del self.all_times
        for output in self.outputs:
            output.raw_output = fastflatten(output.outputarray)
            del output.outputarray
            
    def make_instruction_table(self):
        self.collect_change_times()
        self.make_timeseries()
        self.expand_change_times()
        self.expand_timeseries()
        self.make_raw_output()
 

class Output:
    description = 'generic output'
    # Overridden by subclasses, for example {1:'open', 0:'closed'}
    allowed_states = {}
    
    def __init__(self,name,IO_device,connection_number):
        self.name = name
        self.instructions = {}
        self.connected_to_device = IO_device
        self.connection_number = connection_number
        self.ramp_limits = []
        IO_device.outputs.append(self)
        
        
    def instruction_to_string(self,instruction):
        """gets a human readable description of an instruction"""
        if isinstance(instruction,dict):
            return instruction['description']
        elif self.allowed_states:
            return str(self.allowed_states[instruction])
        else:
            return str(instruction)
        
    def add_instruction(self,time,instruction):
        #Check that this doesn't collide with previous instructions:
        if time in self.instructions.keys():
            err = ' '.join(['WARNING: State of', self.description, self.name, 'at t=%ss'%str(time),
                 'has already been set to %s.'%self.instruction_to_string(self.instructions[time]),
                 'Overwriting to %s.\n'%self.instruction_to_string(instruction)])
            sys.stderr.write(err + '\n')
        # Check that ramps don't collide
        if isinstance(instruction,dict):
            for start, end in self.ramp_limits:
                if start < time < end or start < instruction['end time'] < end:
                    err = ' '.join(['ERROR: State of', self.description, self.name, 'from t = %ss to %ss'%(str(start),str(end)),
                        'has already been set to %s.'%self.instruction_to_string(self.instructions[start]),
                        'Cannot set to %s from t = %ss to %ss. Stopping.'%(self.instruction_to_string(instruction),str(time),str(instruction['end time']))])
                    sys.stderr.write(err + '\n')
                    sys.exit(1)
            self.ramp_limits.append((time,instruction['end time']))
        self.instructions[time] = instruction
        
    def perform_checks(self):
        # Check if there are no instructions. Generate a warning and insert an
        # instruction telling the output to remain at zero.
        if not self.instructions:
            sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions. It will be set to %s for all time.\n'%self.instruction_to_string(0)]))
            self.add_instruction(0,0)  
        # Check if there are no instructions at t=0. Generate a warning and insert an
        # instruction telling the output to start at zero.
        if 0 not in self.instructions.keys():
            sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions at t=0. It will initially be set to %s.\n'%self.instruction_to_string(0)]))
            self.add_instruction(0,0) 
        # Check that ramps have instructions following them.
        # If they don't, insert an instruction telling them to hold their final value.
        for instruction in self.instructions.values():
            if isinstance(instruction, dict) and instruction['end time'] not in self.instructions.keys():
                self.add_instruction(instruction['end time'], instruction['function'](instruction['end time']))
         
    def make_times(self):
        self.perform_checks()
        self.times = self.instructions.keys()
        self.times.sort()
            
    def make_timeseries(self,change_times):
        self.timeseries = []
        i = 0
        for change_time in change_times:
            try:
                if i < len(self.times):
                    while change_time >= self.times[i]:
                        i += 1
            except IndexError:
                # We allow the index to go one higher, since we're
                # intentionally overshooting the mark and are then
                # interested in self.times[i-1].  Raise the error
                # otherwise.
                if not i == len(self.times):
                    raise
            instruction = self.instructions[self.times[i-1]]
            self.timeseries.append(instruction)     
    
    def make_outputarray(self,all_times):
        self.outputarray = []
        for i, time in enumerate(all_times):
            if iterable(time):
                if isinstance(self.timeseries[i],dict):
                    # We evaluate the functions at the midpoints of the
                    # timesteps in order to remove the zero-order hold
                    # error introduced by sampling an analogue signal:
                    try:
                        midpoints = time + 0.5*(time[1] - time[0])
                    except IndexError:
                        # Time array might be only one element long, so we
                        # can't calculate the step size this way. That's
                        # ok, the final midpoint is determined differently
                        # anyway:
                        midpoints = time
                    # We need to know when the first clock tick is after
                    # this ramp ends. It's either an array element or a
                    # single number depending on if this ramp is followed
                    # by another ramp or not:
                    next_time = all_times[i+1][0] if iterable(all_times[i+1]) else all_times[i+1]
                    midpoints[-1] = time[-1] + 0.5*(next_time - time[-1])
                    outarray = self.timeseries[i]['function'](midpoints)
                else:
                    outarray = empty(len(time),dtype=float32)
                    outarray.fill(self.timeseries[i])
                self.outputarray.append(outarray)
            else:
                self.outputarray.append(self.timeseries[i])
    
    def write_raw_output_to_file(self):
        grp = hdf5_file[self.connected_to_device.name]
        grp.create_dataset(str(self.connection_number),data=self.raw_output)
        print 'saved a dataset'
        
         
class NIBoard:
    description = 'NI board'
    clock_limit = 1e6 #TODO get a real number to put here
    def __init__(self, name, clockedby):
        self.name = name
        # Make sure the clock knows is limited by this device's max clock rate:
        if clockedby.clock_limit > self.clock_limit:
            clockedby.clock_limit = self.clock_limit
        self.outputs = clockedby.outputs
        hdf5_file.create_group(self.name)
       
class PulseBlaster(IODevice):
    pb_instructions = {'STOP': 1, 'LOOP': 2, 'END_LOOP': 3}
    description = 'PulseBlaster'
    clock_limit = 25.0e6 # Slight underestimate I think.
    clock_connection = 11
    direct_outputs = []
    
    def perform_checks(self):
        for output in self.outputs:
            if output.connected_to_device is self:
                if not (isinstance(output.connection_number,int) or output.connection_number < 12):
                    sys.stderr.write('%s is set as connected to output connection %d of %s. Output connection \
                                      number must be a integer less than 12\n.'%(output.name, output.connection_number, self.name))
                    sys.exit(1)
                for other_output in self.direct_outputs:
                    if output.connection_number == other_output.connection_number:
                        sys.stderr.write('%s %s and %s %s are both set as connected to output %d of %s! Stopping.\n'%(output.name,
                                         other_output.name, output.connection_number, self.name))
                        sys.exit(1)
                self.direct_outputs.append(output)    
                
    def convert_to_pb_inst(self):
        self.pb_inst = []
        # index to keep track of where in output.raw_output the
        # pulseblaster flags are coming from
        i = 0
        # index to record what line number of the pulseblaster hardware
        # instructions we're up to:
        j = 0
        for instruction in self.clock:
            flags = [0]*12
            for output in self.direct_outputs:
                flags[output.connection_number] = int(output.raw_output[i])
            flags[11] = 1
            flagstring = ''.join([str(flag) for flag in flags])
            if instruction['reps'] > 1048576:
                sys.stderr.write('ERROR: Pulseblaster cannot support more than 1048576 loop iterations. ' +
                                 str(instruction['reps']) +' were requested at t = ' + str(instruction['start']) + '. '+
                                 'This can be fixed easily enough by using nested loops. If it is needed, ' +
                                  'please file a feature request at' +
                                  'http://redmine.physics.monash.edu.au/projects/labscript. Stopping.\n')
                sys.exit(1)
            self.pb_inst.append({'flags': flagstring, 'instruction': 'LOOP',
                                 'data': instruction['reps'], 'delay': instruction['step']*1e9/2.0})
            flags[11] = 0
            flagstring = ''.join([str(flag) for flag in flags])
            self.pb_inst.append({'flags': flagstring, 'instruction': 'END_LOOP',
                                 'data': j, 'delay': instruction['step']*1e9/2.0})
            j += 2
            i += instruction['reps']
        # Gotta put a stop instruction at the end. It will have a short
        # delay time and set everything back to zero:
        self.pb_inst.append({'flags': '000000000000', 'instruction': 'STOP',
                                 'data': 0, 'delay': 10.0/self.clock_limit*1e9})  
                                           
    def write_instructions_to_files(self):
        # The raw output of each analogue output not direcly connected
        # to the PulseBlaster:
        for output in self.outputs:
            if output not in self.direct_outputs:
                output.write_raw_output_to_file()
        # The table of instructions for the PulseBlaster itself:
        pb_dtype = [('freq0',int32), ('phase0',int32), ('amp0',int32), 
                    ('dds_en0',int32), ('phase_reset0',int32),
                    ('freq1',int32), ('phase1',int32), ('amp1',int32),
                    ('dds_en1',int32), ('phase_reset1',int32),
                    ('flags',int32), ('inst',int32),
                    ('inst_data',int32), ('length',float64)]
                    
        pb_inst_table = empty(len(self.pb_inst),dtype = pb_dtype)
        for i,inst in enumerate(self.pb_inst):
            flagint = int(inst['flags'][::-1],2)
            instructionint = self.pb_instructions[inst['instruction']]
            dataint = inst['data']
            delaydouble = inst['delay']
            pb_inst_table[i] = (0,0,0,0,0,0,0,0,0,0, flagint, 
                                instructionint, dataint, delaydouble)
        group = hdf5_file.create_group(self.name)
        group.create_dataset('PULSE_PROGRAM', data = pb_inst_table)
        print 'saved pulse program for %s.'%self.name
        hdf5_file.flush()
        hdf5_file.close()
           
    def generate_code(self):
        self.perform_checks()
        self.make_instruction_table()
        self.convert_to_pb_inst()
        import time
        start_time = time.time()
        self.write_instructions_to_files()
        print "write time:", time.time() - start_time
        


class AnalogueOut(Output):
    description = 'analogue output'
    def ramp(self,t,duration,initial,final,samplerate):
        self.add_instruction(t, {'function': functions.ramp(t,duration,initial,final), 'description':'linear ramp',
                                 'end time': t + duration, 'clock rate': samplerate})
        return duration
    
                                 
    def sine(self,t,duration,amplitude,angfreq,phase,dc_offset,samplerate):
        self.add_instruction(t, {'function': functions.sine(t,amplitude,angfreq,phase,dc_offset), 'description':'sine wave',
                                 'end time': t + duration, 'clock rate': samplerate})
        return duration   
                                 
    def constant(self,t,value):
        self.add_instruction(t,value)
        
        
class DigitalOut(Output):
    description = 'digital output'
    allowed_states = {1:'high', 0:'low'}
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 


class Shutter(DigitalOut):
    description = 'shutter'
    allowed_states = {1:'open', 0:'closed'}  
    def open(self,t):
        self.go_high(t)
    def close(self,t):
        self.go_low(t)


class NovaTechDDS(Output): #make an intermediate device class?
    descpription = 'NovaTech DDS9m'
    def __init__(self):
        pass
    
def DDSOut(Output):
    description='generic DDS device'
    def __init__(self,name,IO_device,connection_number):
        self.name = name
        self.instructions = {}
        self.connected_to_device = IO_device
        self.connection_number = connection_number
        self.ramp_limits = []
        IO_device.outputs.append(self)
    
             
def open_hdf5_file():
    try:
        hdf5_filename = sys.argv[1]
    except:
        sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
        sys.exit(1)
    if not os.path.exists(hdf5_filename):
        sys.stderr.write('ERROR: Provided hdf5 filename %s doesn\'t exist. Stopping.\n'%hdf5_filename)
        sys.exit(1)
    try:
        hdf5_file = h5py.File(hdf5_filename,'a')
    except:
        sys.stderr.write('ERROR: Couldn\'t open %s for writing. '%hdf5_filename +
                         'Check it is a valid hdf5 file and is not read only.\n')
        sys.exit(1) 
    return hdf5_file      
        
                                         
inventory = []
hdf5_file = open_hdf5_file()
params = dict(hdf5_file['params'].attrs)
    
for name in params.keys():
    if name in globals().keys() or name in locals().keys() or name in dir(__builtins__):
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         ' is already a name used by Python, LabScript, or Pylab.'+
                         ' Please choose a different variable name to avoid a conflict.\n')
        sys.exit(1)
    if name in keyword.kwlist:
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         ' is a reserved Python keyword.' +
                         ' Please choose a different variable name.\n')
        sys.exit(1)
        
    try:
        exec(name + ' = 0')
        assert '.' not in name
    except:
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         'is not a valid python variable name.' +
                         ' Please choose a different variable name.\n')
        sys.exit(1)
    exec(name + " = params['%s']"%name )
    















