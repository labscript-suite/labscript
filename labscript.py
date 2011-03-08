import os
import sys
import keyword

import h5py
from pylab import *

import functions

def bitfield(arrays):
    """converts a list of eight arrays of ones and zeros into a single
    array of 8 bit ints"""
    if arrays[0] is 0:
        y = zeros(max([len(arr) if iterable(arr) else 1 for arr in arrays]),dtype=uint8)
    else:
        y = array(arrays[0],dtype=uint8)
    for i in range(1,8):
        y |= arrays[i]<<i
    return y

def fastflatten(inarray, dtype):
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
    flat = empty(total_points,dtype=dtype)
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

def plot_outputs(display=False):
    for device in inventory:
        if device.parent_device is None:
            for i, output in enumerate(device.get_all_outputs()):
                if isinstance(output,Output):
                    t,y = discretise(device.times,output.raw_output/float32(output.scale_factor), device.stop_time)
                    plot(t,y,'-',label=output.name)

    grid(True)
    xlabel('time (seconds)')
    ylabel('output values')
    title('Pseudoclocked outputs')
    legend(loc='upper left')
    #axis([0,max([device.stop_time if isinstance(device, PulseBlaster) else 0 for device in inventory]),-1,5.5])
    if '-show' in sys.argv:
        show()
    else:
        savefig('outputs.png')
        
class Device(object):
    description = 'Generic Device'
    allowed_children = None
    def __init__(self,name,parent_device,connection):
        if self.allowed_children is None:
            allowed_children = [Device]
        self.name = name
        self.parent_device = parent_device
        self.connection = connection
        self.child_devices = []
        inventory.append(self)
        if parent_device:
            parent_device.add_device(self)
        
    def add_device(self,device):
        if any([isinstance(device,DeviceClass) for DeviceClass in self.allowed_children]):
            self.child_devices.append(device)
        else:
            sys.stdout.write('ERROR: devices of type %s cannot be attached to devices of type %s. \n'%(device.description,self.description))
            sys.exit(0)
            
    def get_all_outputs(self):
        all_outputs = []
        for device in self.child_devices:
            if isinstance(device,Output):
                all_outputs.append(device)
            else:
                all_outputs.extend(device.get_all_outputs())
        return all_outputs

    def generate_code(self):
        for device in self.child_devices:
            device.generate_code()
  

class PseudoClock(Device):
    description = 'Generic Pseudoclock'
    allowed_children = [Device]
    def __init__(self,name):
        Device.__init__(self,name,parent_device=None,connection=None)
    
    def collect_change_times(self, outputs):
        """Asks all connected outputs for a list of times that they
        change state. Takes the union of all of these times. Note
        that at this point, a change from holding-a-constant-value
        to ramping-through-values is considered a single state
        change. The clocking times will be filled in later in the
        expand_change_times function, and the ramp values filled in with
        expand_timeseries."""
        change_times = []
        for output in outputs:
            change_times.extend(output.get_change_times())
        # Change to a set and back to get rid of duplicates:
        change_times = list(set(change_times))
        change_times.sort()
        # Check that no two instructions are too close together:
        for i, t in enumerate(change_times[:-1]):
            dt = change_times[i+1] - t
            if dt < 1.0/self.clock_limit:
                sys.stderr.write('ERROR: Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(self.change_times[i+1])) +
                                  'One or more connected devices cannot support update delays shorter than %s sec. Stopping.\n'%str(1.0/self.clock_limit))
                sys.exit(1)
        return change_times
    
    def expand_change_times(self, change_times, outputs):
        """For each time interval delimited by change_times, constructs
        an array of times at which the clock for this device needs to
        tick. If the interval has all outputs having constant values,
        then only the start time is stored.  If one or more outputs are
        ramping, then the clock ticks at the maximum clock rate requested
        by any of the outputs. Also produces a higher level description
        of the clocking; self.clock. This list contains the information
        that facilitates programming a pseudo clock using loops."""
        all_times = []
        clock = []
        for i, time in enumerate(change_times):
            # what's the fastest clock rate?
            maxrate = 0
            for output in outputs:
                # Check if output is sweeping and has highest clock rate
                # so far. If so, store its clock rate to max_rate:
                if isinstance(output.timeseries[i],dict) and output.timeseries[i]['clock rate'] > maxrate:
                    # It does have the highest clock rate? Then store that rate to max_rate:
                    maxrate = output.timeseries[i]['clock rate']
            if maxrate > self.clock_limit:
                sys.stderr.write('ERROR: at t = %s sec, a clock rate of %s Hz was requested. '%(str(time),str(maxrate)) + 
                                 'One or more devices connected to %s cannot support clock rates higher than %sHz. Stopping.\n'%(str(self.name),str(self.clock_limit)))
                sys.exit(1)
                
            if maxrate:
                # If there was ramping at this timestep, how many clock ticks fit before the next instruction?
                n_ticks, remainder = divmod((change_times[i+1] - time)*maxrate,1)
                n_ticks = int(n_ticks)
                # Can we squeeze the final clock cycle in at the end?
                if remainder and remainder/float(maxrate) >= 1/float(self.clock_limit):
                    # Yes we can. Clock speed will be as
                    # requested. Otherwise the final clock cycle will
                    # be too long, by the fraction 'remainder'.
                    n_ticks += 1
                duration = n_ticks/float(maxrate) # avoiding integer division
                ticks = linspace(time,time + duration,n_ticks,endpoint=False)
                all_times.append(array(ticks,dtype=float32))
                # Note that even though all arrays are single precision
                # floating point, the numbers stored to the clock list
                # below are double precision. This is important so that
                # rounding errors in the stepsize don't become significant
                # after many clock cycles.
                if n_ticks > 1:
                    # If n_ticks is only one, then this step doesn't do
                    # anything, it has reps=0. So we should only include
                    # it if n_ticks > 1:
                    clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate)})
                # The last clock tick has a different duration depending on the next step:
                clock.append({'start': ticks[-1], 'reps': 1, 'step': change_times[i+1] - ticks[-1]})
            else:
                all_times.append(time)
                try: 
                    # If there was no ramping, here is a single clock tick:
                    clock.append({'start': time, 'reps': 1, 'step': change_times[i+1] - time})
                except IndexError:
                    if i != len(change_times) - 1:
                        raise
                    # There is no next instruction. Hold the last clock
                    # tick until self.stop_time.
                    if self.stop_time > time:
                        clock.append({'start': time, 'reps': 1, 'step': self.stop_time - time})
                    # Error if self.stop_time has been set to less
                    # than the time of the last instruction:
                    elif self.stop_time < time:
                        sys.stderr.write('ERROR: %s %s has more instructions after the experiment\'s stop time. Stopping.\n'%(self.description,self.name))
                        sys.exit(1)
                    # If self.stop_time is the same as the time of the last
                    # instruction, then we'll get the last instruction
                    # out still, so that the total number of clock
                    # ticks matches the number of data points in the
                    # Output.raw_output arrays. We'll make this last
                    # cycle be at ten times the maximum step duration.
                    else:
                        clock.append({'start': time, 'reps': 1, 'step': 10.0/self.clock_limit})
        return all_times, clock
                        
    def generate_clock(self):
        outputs = self.get_all_outputs()
        change_times = self.collect_change_times(outputs)
        for output in outputs:
            output.make_timeseries(change_times)
        all_times, clock = self.expand_change_times(change_times, outputs)
        for output in outputs:
            output.expand_timeseries(all_times)
        self.clock = clock
        self.times = fastflatten(all_times,float32)
        
    def generate_code(self):
        self.generate_clock()
        Device.generate_code(self)
            

class PulseBlaster(PseudoClock):
    pb_instructions = {'STOP': 1, 'LOOP': 2, 'END_LOOP': 3}
    description = 'PulseBlaster'
    clock_limit = 25.0e6 # Slight underestimate I think.
    clock_connection = 11
    
    def get_direct_outputs(self):
        """Finds out which outputs are directly attached to the PulseBlaster"""
        direct_outputs = []
        for output in self.get_all_outputs():
            if output.parent_device is self:
                if not (isinstance(output.connection,int) and output.connection < 12):
                    sys.stderr.write('%s is set as connected to output connection %d of %s. Output connection \
                                      number must be a integer less than 12\n.'%(output.name, output.connection, self.name))
                    sys.exit(1)
                for other_output in direct_outputs:
                    if output.connection == other_output.connection:
                        sys.stderr.write('%s %s and %s %s are both set as connected to output %d of %s. Stopping.\n'%(output.name,
                                         other_output.name, output.connection, self.name))
                        sys.exit(1)
                direct_outputs.append(output)
        return direct_outputs
    
    def convert_to_pb_inst(self, direct_outputs):
        pb_inst = []
        # index to keep track of where in output.raw_output the
        # pulseblaster flags are coming from
        i = 0
        # index to record what line number of the pulseblaster hardware
        # instructions we're up to:
        j = 0
        for instruction in self.clock:
            flags = [0]*12
            for output in direct_outputs:
                flags[output.connection] = int(output.raw_output[i])
            flags[11] = 1
            flagstring = ''.join([str(flag) for flag in flags])
            if instruction['reps'] > 1048576:
                sys.stderr.write('ERROR: Pulseblaster cannot support more than 1048576 loop iterations. ' +
                                 str(instruction['reps']) +' were requested at t = ' + str(instruction['start']) + '. '+
                                 'This can be fixed easily enough by using nested loops. If it is needed, ' +
                                  'please file a feature request at' +
                                  'http://redmine.physics.monash.edu.au/projects/labscript. Stopping.\n')
                sys.exit(1)
            pb_inst.append({'flags': flagstring, 'instruction': 'LOOP',
                                 'data': instruction['reps'], 'delay': instruction['step']*1e9/2.0})
            flags[11] = 0
            flagstring = ''.join([str(flag) for flag in flags])
            pb_inst.append({'flags': flagstring, 'instruction': 'END_LOOP',
                                 'data': j, 'delay': instruction['step']*1e9/2.0})
            j += 2
            i += instruction['reps']
        # Gotta put a stop instruction at the end. It will have a short
        # delay time and set everything back to zero:
        pb_inst.append({'flags': '000000000000', 'instruction': 'STOP',
                        'data': 0, 'delay': 10.0/self.clock_limit*1e9})  
        # OK now we squeeze the instructions into a numpy array ready for writing to hdf5:
        pb_dtype = [('freq0',int32), ('phase0',int32), ('amp0',int32), 
                    ('dds_en0',int32), ('phase_reset0',int32),
                    ('freq1',int32), ('phase1',int32), ('amp1',int32),
                    ('dds_en1',int32), ('phase_reset1',int32),
                    ('flags',int32), ('inst',int32),
                    ('inst_data',int32), ('length',float64)]
        pb_inst_table = empty(len(pb_inst),dtype = pb_dtype)
        for i,inst in enumerate(pb_inst):
            flagint = int(inst['flags'][::-1],2)
            instructionint = self.pb_instructions[inst['instruction']]
            dataint = inst['data']
            delaydouble = inst['delay']
            pb_inst_table[i] = (0,0,0,0,0,0,0,0,0,0, flagint, 
                                instructionint, dataint, delaydouble)
                              
        # Okey now write it to the file:   
        group = hdf5_file.create_group(self.name)
        group.create_dataset('PULSE_PROGRAM', compression=compression,data = pb_inst_table)         
                              
    def generate_code(self):
        PseudoClock.generate_code(self)
        direct_outputs = self.get_direct_outputs()
        self.convert_to_pb_inst(direct_outputs)


            
class Output(Device):
    description = 'generic output'
    allowed_states = {}
    dtype = float32
    scale_factor = 1
    def __init__(self,name,parent_device,connection):
        self.instructions = {}
        self.ramp_limits = [] # For checking ramps don't overlap
        Device.__init__(self,name,parent_device,connection)      

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
        
    def get_change_times(self):
        """If this function is being called, it means that the parent
        Pseudoclock has requested a list of times that this output changes
        state. First we'll need to perform some checks to make sure that
        the instructions the user has entered make sense. Then the list
        of times is returned."""
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
        times = self.instructions.keys()
        self.times = times
        return times
        
    def make_timeseries(self,change_times):
        """If this is being called, then it means the parent Pseudoclock
        has asked for a list of this output's states at each time in
        change_times. (Which are the times that one or more connected
        outputs in the same pseudoclock change state). By state, I don't
        mean the value of the output at that moment, rather I mean what
        instruction it has. This might be a single value, or it might
        be a reference to a function for a ramp etc. This list of states
        is stored in self.timeseries rather than being returned."""
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
        
    def expand_timeseries(self,all_times):
        """This function evaluates the ramp functions in self.timeseries
        at the time points in all_times, and creates an array of output
        values at those times.  These are the values that this output
        should update to on each clock tick, and are the raw values that
        should be used to program the output device.  They are stored
        in self.raw_output."""
        outputarray = []
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
                outputarray.append(outarray)
            else:
                outputarray.append(self.timeseries[i])
        del self.timeseries # don't need this any more.
        self.raw_output = fastflatten(outputarray, self.dtype)
        

class AnalogueOut(Output):
    description = 'analogue output'
    def ramp(self,t,duration,initial,final,samplerate):
        self.add_instruction(t, {'function': functions.ramp(t,duration,initial,final), 'description':'linear ramp',
                                 'end time': t + duration, 'clock rate': samplerate})
        return duration
    
                                 
    def sine(self,t,duration,amplitude,angfreq,phase,dc_offset,samplerate):
        self.add_instruction(t, {'function': functions.sine(t,duration,amplitude,angfreq,phase,dc_offset), 'description':'sine wave',
                                 'end time': t + duration, 'clock rate': samplerate})
        return duration   
                                 
    def constant(self,t,value):
        self.add_instruction(t,value)
        
        
class DigitalOut(Output):
    description = 'digital output'
    allowed_states = {1:'high', 0:'low'}
    dtype = uint8
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 


class NIBoard(Device):
    allowed_children = [AnalogueOut, DigitalOut]
    n_analogues = 4
    n_digiports = 2 # number of 'ports', 8 digital outputs per port.
    clock_limit = 500e3 # underestimate I think.
    
    def __init__(self, name, parent_device):
        Device.__init__(self, name, parent_device, connection=None)
        self.parent_device.clock_limit = min([self.parent_device.clock_limit,self.clock_limit])
        
    def add_output(self,output):
        # TODO: check there are no duplicates, check that connection
        # string is formatted correctly.
        Device.add_output(self,output)
        
    def convert_bools_to_bytes(self,digitals):
        ports = {}
        bitfields = {}
        for i in range(self.n_digiports):
            ports['p%d'%i] = [0]*8
        for output in digitals:
            port, line = output.connection.strip('p').split('.')
            port, line  = int(port),int(line)
            ports['p%d'%port][line] = output.raw_output
        for port, data in ports.items():
            bitfields[port] = bitfield(data)
        return bitfields
            
            
    def convert_to_int16(self,output):
        """converts floats between -10 and 10 to signed 16 bit integers
        between -32767 and 32767 inclusive (excluding -32768 for symmetry)"""
        # This function is currenly responsible for about a 50% increase
        # in run time, but i'm not sure whether it can be optimised much
        # without sacrificing correctness.  The any(output.raw_output >
        # 10) checks etc are what are slow, but without bounds checking
        # the type casting to int16 wraps around without detecting
        # overflow. Writing int16s results in less time to write to
        # disk though, and so the net effect is not much. However if
        # we use an OS with proper caching, such as Win 7, then write
        # times don't concern us much. If we pass float32's, then LabVIEW
        # presumably does its own bounds checking before giving the data
        # to its cards. If it doesn't then we'll have to do this anyway,
        # grumble grumble, to make sure the experiment doesn't silently
        # fail when the user enters too high a voltage...
        if any(output.raw_output > 10 )  or any(output.raw_output < -10 ):
            sys.stderr.write('ERROR: %s %s '%(output.description, output.name) +
                              'can only have values between -10 and 10 Volts, ' + 
                              'the limit imposed by %s. Stopping.\n'%self.name)
            sys.exit(1)
        # have to round first else we just get the floor function on absolute values:
        output.raw_output = array((3276.7*output.raw_output).round(),dtype=int16)
        output.scale_factor = 3276.7
    
    def generate_code(self):
        analogues = {}
        digitals = {}
        for output in self.child_devices:
            if isinstance(output,AnalogueOut):
                analogues[output.connection] = output
            else:
                digitals[output.connection] = output
        digiports = self.convert_bools_to_bytes(digitals.values())
        dtypes = [('ao%d'%i,int16) for i in range(self.n_analogues)] + \
                 [('p%d'%i,uint8) for i in range(self.n_digiports)]
        out_table = zeros(len(self.parent_device.times),dtype=dtypes)
        for connection, analogueout in analogues.items():
            self.convert_to_int16(analogueout)
            out_table[connection] = analogueout.raw_output
        for portno, data in digiports.items():
            out_table[portno] = data
        grp = hdf5_file.create_group(self.name)
        grp.attrs['analogue scale factor'] = 3276.7
        grp.create_dataset('OUTPUT_VALUES',compression=compression,data=out_table)


class Shutter(DigitalOut):
    description = 'shutter'
    allowed_states = {1:'open', 0:'closed'}  
    def open(self,t):
        self.go_high(t)
    def close(self,t):
        self.go_low(t)  
  
class DDS(Device):
    description = 'DDS output'
    allowed_children = [AnalogueOut] # Adds its own children when initialised
    def __init__(self,name,parent_device,connection):
        Device.__init__(self,name,parent_device,connection)
        self.frequency = AnalogueOut(self.name+' (freq)',self,None)
        self.amplitude = AnalogueOut(self.name+' (amp)',self,None)
        self.phase = AnalogueOut(self.name+' (phase)',self,None)
    def setamp(self,t,value):
        self.amplitude.constant(t,value)
    def setfreq(self,t,value):
        self.frequency.constant(t,value)
    def setphase(self,t,value):
        self.phase.constant(t,value)
        
class NovaTechDDS9M(Device):
    description = 'Novatech DDS'
    allowed_children = [DDS]
    clock_limit = 500e3 # TODO: find out what the actual max clock rate is.
    
    def __init__(self,name,parent_device):
        Device.__init__(self,name,parent_device,None)
        self.parent_device.clock_limit = min([self.parent_device.clock_limit,self.clock_limit])
        
    def quantise_freq(self,output):
        # Ensure that frequencies are within bounds:
        if any(output.raw_output > 171e6 )  or any(output.raw_output < 0.1 ):
            sys.stderr.write('ERROR: %s %s '%(output.description, output.name) +
                              'can only have values between 0.1Hz and 171MHz, ' + 
                              'the limit imposed by %s. Stopping.\n'%self.name)
            sys.exit(1)
        # It's faster to add 0.5 then typecast than to round to integers first:
        output.raw_output = array((1e7*output.raw_output)+0.5,dtype=uint32)
        output.scale_factor = 1e7
        
    def quantise_phase(self,output):
        # ensure that phase wraps around:
        output.raw_output %= 360
        # It's faster to add 0.5 then typecast than to round to integers first:
        output.raw_output = array((45.511111111111113*output.raw_output)+0.5,dtype=uint16)
        output.scale_factor = 45.511111111111113
        
    def quantise_amp(self,output):
        # ensure that amplitudes are within bounds:
        if any(output.raw_output > 1 )  or any(output.raw_output < 0):
            sys.stderr.write('ERROR: %s %s '%(output.description, output.name) +
                              'can only have values between 0 and 1 (Volts peak to peak approx), ' + 
                              'the limit imposed by %s. Stopping.\n'%self.name)
        # It's faster to add 0.5 then typecast than to round to integers first:
        output.raw_output = array((1023*output.raw_output)+0.5,dtype=uint16)
        output.scale_factor = 1023
        
    def generate_code(self):
        DDSs = {}
        for output in self.child_devices:
        # Check that the instructions will fit into RAM:
            if len(output.frequency.raw_output) > 32768:
                sys.stderr.write('ERROR: %s can only support 32768 instructions. '%self.name +
                                 'Please decrease the sample rates of devices on the same clock, ' + 
                                 'or connect %s to a different pseudoclock. Stopping.\n'%self.name)
                sys.exit(1)
            # TODO: check that there are only two and that their connection no.s are right.
            DDSs[output.connection] = output
        for dds in DDSs.values():
            self.quantise_freq(dds.frequency)
            self.quantise_phase(dds.phase)
            self.quantise_amp(dds.amplitude)
            
        dtypes = [('freq%d'%i,uint32) for i in range(2)] + \
                 [('phase%d'%i,uint16) for i in range(2)] + \
                 [('amp%d'%i,uint16) for i in range(2)]
        out_table = zeros(len(self.parent_device.times),dtype=dtypes)
        for connection, dds in DDSs.items():
            out_table['freq%d'%connection] = dds.frequency.raw_output
            out_table['amp%d'%connection] = dds.amplitude.raw_output
            out_table['phase%d'%connection] = dds.phase.raw_output
        grp = hdf5_file.create_group(self.name)
        grp.attrs['frequency scale factor'] = 1e7
        grp.attrs['amplitude scale factor'] = 1023
        grp.attrs['phase scale factor'] = 45.511111111111113
        grp.create_dataset('TABLE_DATA',compression=compression,data=out_table) 
        
def stop(t):
    for device in inventory:
        if not device.parent_device:  
            device.stop_time = t
            
def generate_code():
    for device in inventory:
        if not device.parent_device:
            device.generate_code()
            print
            print device.name, '\t', len(device.times), 'x', device.times.dtype 
            for output in device.get_all_outputs():
                print output.name, '\t', len(output.raw_output), 'x', output.raw_output.dtype
            print
    hdf5_file.close()
 
def open_hdf5_file():
    try:
        hdf5_filename = sys.argv[-1]
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
if '-compress' in sys.argv:
    compression = 'gzip'
else:
    compression = None
        
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
    


       
