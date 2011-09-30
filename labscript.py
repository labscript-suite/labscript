import __main__

import os
import sys
import keyword
import gtk

import h5py
from pylab import *

import functions

def bitfield(arrays,dtype):
    """converts a list of arrays of ones and zeros into a single
    array of unsigned ints of the given datatype."""
    n = {uint8:8,uint16:16,uint32:32}
    if arrays[0] is 0:
        y = zeros(max([len(arr) if iterable(arr) else 1 for arr in arrays]),dtype=dtype)
    else:
        y = array(arrays[0],dtype=dtype)
    for i in range(1,n[dtype]):
        if iterable(arrays[i]):
            print arrays[i], arrays[i]<<i
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
        try:
            # Test that name is a valid Python variable name:
            exec '%s = None'%name
            assert '.' not in name
        except:
            raise ValueError('%s is not a valid Python variable name. Stopping.\n'%name)
        # Put self into the global namespace:
        setattr(__main__,name,self)
        
    def add_device(self,device):
        if any([isinstance(device,DeviceClass) for DeviceClass in self.allowed_children]):
            self.child_devices.append(device)
        else:
            sys.stdout.write('ERROR: devices of type %s cannot be attached to devices of type %s. \n'%(device.description,self.description))
            sys.exit(1)
            
    def get_all_outputs(self):
        all_outputs = []
        for device in self.child_devices:
            if isinstance(device,Output):
                all_outputs.append(device)
            else:
                all_outputs.extend(device.get_all_outputs())
        return all_outputs
    
    def get_all_children(self):
        all_children = []
        for device in self.child_devices:
              all_children.append(device)
              all_children.extend(device.get_all_children())
        return all_children

    def generate_code(self):
        for device in self.child_devices:
            device.generate_code()
  

class PseudoClock(Device):
    description = 'Generic Pseudoclock'
    allowed_children = [Device]
    generation = 0
    def __init__(self,name,parent_device=None,connection=None):
        Device.__init__(self,name,parent_device,connection)
    
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
                sys.stderr.write('ERROR: Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(change_times[i+1])) +
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
                    # it if n_ticks > 1.
                    if n_ticks > 2:
                        #If there is more than one clock tick here,
                        #then we split the ramp into an initial clock
                        #tick, during which the slow clock ticks, and
                        #the rest of the ramping time, during which the
                        #slow clock does not tick.
                        clock.append({'start': time, 'reps': 1, 'step': 1/float(maxrate),'slow_clock_tick':True})
                        clock.append({'start': time + 1/float(maxrate), 'reps': n_ticks-2, 'step': 1/float(maxrate),'slow_clock_tick':False})
                    else:
                        clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate),'slow_clock_tick':True})
                # The last clock tick has a different duration depending
                # on the next step. The slow clock must tick here if it
                # hasn't ticked already, that is if n_ticks = 1.
                clock.append({'start': ticks[-1], 'reps': 1, 'step': change_times[i+1] - ticks[-1],'slow_clock_tick': True if n_ticks == 1 else False})
            else:
                all_times.append(time)
                try: 
                    # If there was no ramping, here is a single clock tick:
                    clock.append({'start': time, 'reps': 1, 'step': change_times[i+1] - time,'slow_clock_tick':True})
                except IndexError:
                    if i != len(change_times) - 1:
                        raise
                    # There is no next instruction. Hold the last clock
                    # tick until self.stop_time.
                    if self.stop_time > time:
                        clock.append({'start': time, 'reps': 1, 'step': self.stop_time - time,'slow_clock_tick':True})
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
                        clock.append({'start': time, 'reps': 1, 'step': 10.0/self.clock_limit,'slow_clock_tick':True})
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
        self.change_times = fastflatten(change_times, float32)
        self.times = fastflatten(all_times,float32)
        
    def generate_code(self):
        self.generate_clock()
        Device.generate_code(self)
            

class PulseBlaster(PseudoClock):
    pb_instructions = {'CONTINUE':0,
                       'STOP': 1, 
                       'LOOP': 2, 
                       'END_LOOP': 3,
                       'BRANCH': 6,
                       'WAIT': 8}
    description = 'PB-DDSII-300'
    clock_limit = 8.3e6 # Slight underestimate I think.
    fast_clock_flag = 0
    slow_clock_flag = 1
    clock_type = 'slow clock'
    
    def get_direct_outputs(self):
        """Finds out which outputs are directly attached to the PulseBlaster"""
        dig_outputs = []
        dds_outputs = []
        for output in self.get_all_outputs():
            # If the device's parent is a DDS (remembering that DDSs
            # have three fake child devices for amp, freq and phase),
            # then maybe that DDS is one of our direct outputs:
            if isinstance(output.parent_device,DDS) and output.parent_device.parent_device is self:
                # If this is the case, then we're interested in that DDS. But we don't want to count it three times:
                if not output.parent_device in dds_outputs:
                    output = output.parent_device
            if output.parent_device is self:
                try:
                    prefix, connection = output.connection.split()
                    assert prefix == 'flag' or prefix == 'dds'
                    connection = int(connection)
                except:
                    sys.stderr.write('%s %s has invalid connection string: \'%s\'. '%(output.description,output.name,str(output.connection)) + 
                                     'Format must be \'flag n\' with n an integer less than 12, or \'dds n\' with n less than 2. Stopping.\n')
                    sys.exit(1)
                if not connection < 12:
                    sys.stderr.write('%s is set as connected to output connection %d of %s. Output connection \
                                      number must be a integer less than 12\n.'%(output.name, connection, self.name))
                    sys.exit(1)
                if prefix == 'dds' and not connection < 2:
                    sys.stderr.write('%s is set as connected to output connection %d of %s. DDS output connection \
                                      number must be a integer less than 2\n.'%(output.name, connection, self.name))
                    sys.exit(1)
                if prefix == 'flag' and connection in [self.slow_clock_flag, self.fast_clock_flag]:
                    sys.stderr.write('%s is set as connected to flag %d of %s.'%(output.name, connection, self.name) +
                                     'This is one of the PulseBlaster\'s clock flags. Stopping.\n')
                    sys.exit(1)
                for other_output in dig_outputs + dds_outputs:
                    if output.connection == other_output.connection:
                        sys.stderr.write('%s %s and %s %s are both set as connected to %s of %s. Stopping.\n'%(output.name,
                                         other_output.name, output.connection, self.name))
                        sys.exit(1)
                if isinstance(output,DigitalOut):
                	dig_outputs.append(output)
                elif isinstance(output, DDS):
                	dds_outputs.append(output)
                
        return dig_outputs, dds_outputs

    def generate_registers(self,dds_outputs):
        ampdicts = {}
        phasedicts = {}
        freqdicts = {}
        group = hdf5_file.create_group('/devices/'+self.name)
        for output in dds_outputs:
            num = int(output.connection.split()[1])
            
            # Ensure that amplitudes are within bounds:
            if any(output.amplitude.raw_output > 1)  or any(output.amplitude.raw_output < 0):
                sys.stderr.write('ERROR: %s %s '%(output.amplitude.description, output.amplitude.name) +
                                  'can only have values between 0 and 1, ' + 
                                  'the limit imposed by %s. Stopping.\n'%output.name)
                                  
            # Ensure that frequencies are within bounds:
            if any(output.frequency.raw_output > 100e6 )  or any(output.frequency.raw_output < 0):
                sys.stderr.write('ERROR: %s %s '%(output.frequency.description, output.frequency.name) +
                                  'can only have values between 0Hz and and 100MHz, ' + 
                                  'the limit imposed by %s. Stopping.\n'%output.name)
                                  
            # Ensure that phase wraps around:
            output.phase.raw_output %= 360
            
            amps = set(output.amplitude.raw_output)
            phases = set(output.phase.raw_output)
            freqs = set(output.frequency.raw_output)
                                  
            if len(amps) > 1024:
                sys.stderr.write('%s dds%d can only support 1024 amplitude registers, and %s have been requested. Stopping.\n'%(self.name, num, str(len(amps))))
                sys.exit(1)
            if len(phases) > 128:
                sys.stderr.write('%s dds%d can only support 128 phase registers, and %s have been requested. Stopping.\n'%(self.name, num, str(len(phases))))
                sys.exit(1)
            if len(freqs) > 1024:
                sys.stderr.write('%s dds%d can only support 1024 frequency registers, and %s have been requested. Stopping.\n'%(self.name, num, str(len(freqs))))
                sys.exit(1)
                
            # start counting at 1 to leave room for the dummy instruction,
            # which LabVEIW will fill in with the state of the front
            # panel:
            ampregs = range(1,len(amps)+1)
            freqregs = range(1,len(freqs)+1)
            phaseregs = range(1,len(phases)+1)
            
            ampdicts[num] = dict(zip(amps,ampregs))
            freqdicts[num] = dict(zip(freqs,freqregs))
            phasedicts[num] = dict(zip(phases,phaseregs))
            
            # The zeros are the dummy instructions:
            freq_table = array([0] + sorted(freqs), dtype = float64)
            amp_table = array([0] + sorted(amps), dtype = float32)
            phase_table = array([0] + sorted(phases), dtype = float64)
            
            subgroup = group.create_group('DDS%d'%num)
            subgroup.create_dataset('FREQ_REGS', compression=compression,data = freq_table)
            subgroup.create_dataset('AMP_REGS', compression=compression,data = amp_table)
            subgroup.create_dataset('PHASE_REGS', compression=compression,data = phase_table)
            
        return freqdicts, ampdicts, phasedicts
        
    def convert_to_pb_inst(self, dig_outputs, dds_outputs, freqs, amps, phases):
        pb_inst = []
        # index to keep track of where in output.raw_output the
        # pulseblaster flags are coming from
        i = 0
        # index to record what line number of the pulseblaster hardware
        # instructions we're up to:
        j = 0
        # We've delegated the initial instruction off to LabVIEW, which
        # can ensure continuity with the state of the front panel. Thus
        # this first instruction doesn't actually do anything.  An initial
        # instruction with the fast clock at zero and the slow clock
        # at one:
        flags = [0]*12
        freqregs = [0]*2
        ampregs = [0]*2
        phaseregs = [0]*2
        flags[self.fast_clock_flag] = 0
        flags[self.slow_clock_flag] = 1 
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs,
                        'flags': ''.join([str(flag) for flag in flags]), 'instruction': 'STOP',
                        'data': 0, 'delay': 10.0/self.clock_limit*1e9})  
        j += 1
        for k, instruction in enumerate(self.clock):
            flags = [0]*12
            freqregs = [0]*2
            ampregs = [0]*2
            phaseregs = [0]*2
            for output in dig_outputs:
                flagindex = int(output.connection.split()[1])
                flags[flagindex] = int(output.raw_output[i])
            for output in dds_outputs:
                ddsnumber = int(output.connection.split()[1])
                freqregs[ddsnumber] = freqs[ddsnumber][output.frequency.raw_output[i]]
                ampregs[ddsnumber] = amps[ddsnumber][output.amplitude.raw_output[i]]
                phaseregs[ddsnumber] = phases[ddsnumber][output.phase.raw_output[i]]
                
            flags[self.fast_clock_flag] = 1
            flags[self.slow_clock_flag] = 0 if instruction['slow_clock_tick'] else 1
            flagstring = ''.join([str(flag) for flag in flags])
            if instruction['reps'] > 1048576:
                sys.stderr.write('ERROR: Pulseblaster cannot support more than 1048576 loop iterations. ' +
                                 str(instruction['reps']) +' were requested at t = ' + str(instruction['start']) + '. '+
                                 'This can be fixed easily enough by using nested loops. If it is needed, ' +
                                  'please file a feature request at' +
                                  'http://redmine.physics.monash.edu.au/projects/labscript. Stopping.\n')
                sys.exit(1)
            pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs,
                            'flags': flagstring, 'instruction': 'LOOP',
                            'data': instruction['reps'], 'delay': instruction['step']*1e9/2.0})
            flags[self.fast_clock_flag] = 0
            flags[self.slow_clock_flag] = 1
            flagstring = ''.join([str(flag) for flag in flags])
            pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs,
                            'flags': flagstring, 'instruction': 'END_LOOP',
                            'data': j, 'delay': instruction['step']*1e9/2.0})
            j += 2
            try:
                if self.clock[k+1]['slow_clock_tick']:
                    i += 1
            except IndexError:
                pass
        # This is how we stop the pulse program. We wait on the second
        # last instruction whilst the control software in LabVIEW
        # reprograms the zeroth instruction to whatever the system
        # is going to revert to once the experiment is done. Then a
        # software trigger is given to the PulseBlaster telling it
        # to resume execution. The program proceeds to the last instruction,
        # which is a branch to zero, and the system ends up in the state
        # of LabVIEW's front panel without missing a beat.
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs,
                        'flags': flagstring, 'instruction': 'WAIT',
                        'data': 0, 'delay': 10.0/self.clock_limit*1e9})  
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs,
                        'flags': flagstring, 'instruction': 'BRANCH',
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
            freq0 = inst['freqs'][0]
            freq1 = inst['freqs'][1]
            phase0 = inst['phases'][0]
            phase1 = inst['phases'][1]
            amp0 = inst['amps'][0]
            amp1 = inst['amps'][1]
            pb_inst_table[i] = (freq0,phase0,amp0,1,0,freq1,amp1,phase1,1,0, flagint, 
                                instructionint, dataint, delaydouble)
                              
        # Okey now write it to the file: 
        group = hdf5_file['/devices/'+self.name]  
        group.create_dataset('PULSE_PROGRAM', compression=compression,data = pb_inst_table)         
        group.create_dataset('FAST_CLOCK', compression=compression,data = self.times)         
        group.create_dataset('SLOW_CLOCK', compression=compression,data = self.change_times)         
#        for thing in pb_inst:
#            for key,val in thing.items():
#                print str(val).center(15),
#            print
                              
    def generate_code(self):
        PseudoClock.generate_code(self)
        dig_outputs, dds_outputs = self.get_direct_outputs()
        freqs, amps, phases = self.generate_registers(dds_outputs)
        self.convert_to_pb_inst(dig_outputs, dds_outputs, freqs, amps, phases)


            
class Output(Device):
    description = 'generic output'
    allowed_states = {}
    dtype = float32
    scale_factor = 1
    generation = 3
    def __init__(self,name,parent_device,connection):
        self.instructions = {}
        self.ramp_limits = [] # For checking ramps don't overlap
        self.clock_type = parent_device.clock_type
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
        #Check that time is not negative:
        if time < 0:
            err = ' '.join(['ERROR:', self.description, self.name, 'has an instruction at t=%ss.'%str(time),
                 'Labscript cannot handle instructions earlier than the experiment start time (t=0). Stopping.\n'])
            sys.stderr.write(err + '\n')
            sys.exit(1)
        #Check that this doesn't collide with previous instructions:
        if time in self.instructions.keys():
            err = ' '.join(['WARNING: State of', self.description, self.name, 'at t=%ss'%str(time),
                 'has already been set to %s.'%self.instruction_to_string(self.instructions[time]),
                 'Overwriting to %s.\n'%self.instruction_to_string(instruction)])
            sys.stderr.write(err + '\n')
        # Check that ramps don't collide
        if isinstance(instruction,dict):
            # No ramps allowed if this output is on a slow clock:
            if self.clock_type == 'slow clock':
                sys.stderr.write('ERROR: %s %s is on a slow clock.'%(self.description, self.name) + 
                                 'It cannot have a function ramp as an instruction. Stopping.\n')
                sys.exit(1)
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
            sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions. It will be set to %s for all time.\n'%self.instruction_to_string(self.default_value)]))
            self.add_instruction(0,self.default_value)  
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
        times.sort()
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
        # If this output is on the slow clock, then its timeseries should
        # not be expanded. It's already as expanded as it'll get.
        if self.clock_type == 'slow clock':
            self.raw_output = fastflatten(self.timeseries,self.dtype)
            return
        outputarray = []
        for i, time in enumerate(all_times):
            if iterable(time):
                if isinstance(self.timeseries[i],dict):
                    # We evaluate the functions at the midpoints of the
                    # timesteps in order to remove the zero-order hold
                    # error introduced by sampling an analog signal:
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
        

class AnalogOut(Output):
    description = 'analog output'
    default_value = 0
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
    default_value = 0
    dtype = uint32
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 


class AnalogIn(Device):
    description = 'Analog Input'
    generation = 3
    def __init__(self,name,parent_device,connection,scale_factor=1.0,units='Volts'):
         self.acquisitions = []
         self.scale_factor = scale_factor
         self.units=units
         Device.__init__(self,name,parent_device,connection)
   
    def acquire(self,label,start_time,end_time,scale_factor=None,units=None):
        if scale_factor is None:
            scale_factor = self.scale_factor
        if units is None:
            units = self.units
        self.acquisitions.append({'start_time': start_time, 'end_time': end_time,
                                 'label': label, 'scale_factor':scale_factor,'units':units})
        
  
class IntermediateDevice(Device):
    generation = 1
    def __init__(self, name, parent_device,clock_type):
        if not clock_type in ['fast clock', 'slow clock']:
            sys.stderr.write('Clock type for %s %s can only be \'slow clock\' or \'fast clock\'. Stopping\n'%(self.name,self.description))
            sys.exit(1)
        Device.__init__(self, name, parent_device, clock_type)
        self.clock_type = clock_type
        self.parent_device.clock_limit = min([self.parent_device.clock_limit,self.clock_limit])
        
        
class NIBoard(IntermediateDevice):
    allowed_children = [AnalogOut, DigitalOut, AnalogIn]
    n_analogs = 4
    n_digitals = 32
    digital_dtype = uint32
    clock_limit = 500e3 # underestimate I think.
    description = 'generic_NI_Board'
    def __init__(self, name, parent_device,clock_type,acquisition_rate=0):
        IntermediateDevice.__init__(self, name, parent_device,clock_type)
        self.acquisition_rate = acquisition_rate
        
    def add_output(self,output):
        # TODO: check there are no duplicates, check that connection
        # string is formatted correctly.
        Device.add_output(self,output)
        
    def convert_bools_to_bytes(self,digitals):
        """converts digital outputs to an array of bitfields stored
        as self.digital_dtype"""
        outputarray = [0]*self.n_digitals
        for output in digitals:
            port, line = output.connection.replace('port','').replace('line','').split('/')
            port, line  = int(port),int(line)
            if port > 0:
                sys.stderr.write('ERROR: Ports > 0 on NI Boards not implemented. Please use port 0, or file a feature request at redmine.physics.monash.edu.au/labscript. Stopping.\n')
                sys.exit(1)
            outputarray[line] = output.raw_output
        bits = bitfield(outputarray,dtype=self.digital_dtype)
        return bits
            
            
    def generate_code(self):
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
        analog_out_table = empty((len(self.parent_device.times),len(analogs)), dtype=float32)
        analog_connections = analogs.keys()
        analog_connections.sort()
        analog_out_attrs = []
        for i, connection in enumerate(analog_connections):
            output = analogs[connection]
            if any(output.raw_output > 10 )  or any(output.raw_output < -10 ):
                # Bounds checking:
                sys.stderr.write('ERROR: %s %s '%(output.description, output.name) +
                                  'can only have values between -10 and 10 Volts, ' + 
                                  'the limit imposed by %s. Stopping.\n'%self.name)
                sys.exit(1)
            analog_out_table[:,i] = output.raw_output
            analog_out_attrs.append(self.name+'/'+connection)
        input_connections = inputs.keys()
        input_connections.sort()
        input_attrs = []
        acquisitions = []
        for connection in input_connections:
            input_attrs.append(self.name+'/'+connection)
            for acq in inputs[connection].acquisitions:
                acquisitions.append((connection,acq['label'],acq['start_time'],acq['end_time'],acq['scale_factor'],acq['units']))
        # The 'a256' dtype below limits the string fields to 256
        # characters. Can't imagine this would be an issue, but to not
        # specify the string length (using dtype=str) causes the strings
        # to all come out empty.
        acquisitions_table_dtypes = [('connection','a256'), ('label','a256'), ('start',float),
                                     ('stop',float), ('scale factor',float), ('units','a256')]
        acquisition_table= empty(len(acquisitions), dtype=acquisitions_table_dtypes)
        for i, acq in enumerate(acquisitions):
            acquisition_table[i] = acq

        digital_out_table = self.convert_bools_to_bytes(digitals.values())
        
        grp = hdf5_file.create_group('/devices/'+self.name)
        analog_dataset = grp.create_dataset('ANALOG_OUTS',compression=compression,data=analog_out_table)
        digital_dataset = grp.create_dataset('DIGITAL_OUTS',compression=compression,data=digital_out_table)
        if len(acquisition_table) > 0:
            input_dataset = grp.create_dataset('ACQUISITIONS',compression=compression,data=acquisition_table)
        grp.attrs['analog_out_channels'] = ', '.join(analog_out_attrs)
        grp.attrs['analog_in_channels'] = ', '.join(input_attrs)
        grp.attrs['digital_lines'] = '/'.join((self.name,'port0','line0:%d'%(self.n_digitals-1)))
        grp.attrs['acquisition_rate'] = self.acquisition_rate
        
class NI_PCI_6733(NIBoard):
    description = 'NI-PCI-6733'
    n_analogs = 8
    n_digitals = 8
    digital_dtype = uint8
    
class NI_PCIe_6363(NIBoard):
    description = 'NI-PCIe-6363'
    n_analogs = 4
    n_digitals = 32
    digital_dtype = uint32
    
class Shutter(DigitalOut):
    description = 'shutter'
    allowed_states = {1:'open', 0:'closed'}  
    def __init__(self,name,parent_device,connection,delay=(0,0)):
        DigitalOut.__init__(self,name,parent_device,connection)
        self.open_delay, self.close_delay = delay
        
    # If a shutter is asked to do something at t=0, it cannot start moving
    # earlier than that.  So initial shutter states will have imprecise
    # timing. Not throwing a warning here because if I did, every run
    # would throw a warning for every shutter. The documentation will
    # have to make a point of this.
    def open(self,t):
        self.go_high(t-self.open_delay if t >= self.open_delay else 0)
    def close(self,t):
        self.go_low(t-self.close_delay if t >= self.close_delay else 0)  
  
class DDS(Device):
    description = 'DDS'
    allowed_children = [AnalogOut] # Adds its own children when initialised
    generation = 2
    def __init__(self,name,parent_device,connection):
        self.clock_type = parent_device.clock_type
        Device.__init__(self,name,parent_device,connection)
        self.frequency = AnalogOut(self.name+'_freq',self,None)
        self.amplitude = AnalogOut(self.name+'_amp',self,None)
        self.phase = AnalogOut(self.name+'_phase',self,None)
        if isinstance(self.parent_device,NovaTechDDS9M):
            self.frequency.default_value = 0.1
    def setamp(self,t,value):
        self.amplitude.constant(t,value)
    def setfreq(self,t,value):
        self.frequency.constant(t,value)
    def setphase(self,t,value):
        self.phase.constant(t,value)
        
class NovaTechDDS9M(IntermediateDevice):
    description = 'NT-DDS9M'
    allowed_children = [DDS]
    clock_limit = 500e3 # TODO: find out what the actual max clock rate is.
    
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
            sys.exit(1)
        # It's faster to add 0.5 then typecast than to round to integers first:
        output.raw_output = array((1023*output.raw_output)+0.5,dtype=uint16)
        output.scale_factor = 1023
        
    def generate_code(self):
        DDSs = {}
        for output in self.child_devices:
        # Check that the instructions will fit into RAM:
            if len(output.frequency.raw_output) > 16383:
                sys.stderr.write('ERROR: %s can only support 16383 instructions. '%self.name +
                                 'Please decrease the sample rates of devices on the same clock, ' + 
                                 'or connect %s to a different pseudoclock. Stopping.\n'%self.name)
                sys.exit(1)
            # TODO: check that there are only two and that their connection no.s are right.
            try:
                prefix, channel = output.connection.split()
                channel = int(channel)
            except:
                sys.stderr.write('%s %s has invalid connection string: \'%s\'. '%(output.description,output.name,str(output.connection)) + 
                                     'Format must be \'channel n\' with n either 0 or 1. Stopping.\n')
                sys.exit(1)
            DDSs[channel] = output
        for dds in DDSs.values():
            self.quantise_freq(dds.frequency)
            self.quantise_phase(dds.phase)
            self.quantise_amp(dds.amplitude)
            
        dtypes = [('freq%d'%i,uint32) for i in range(2)] + \
                 [('phase%d'%i,uint16) for i in range(2)] + \
                 [('amp%d'%i,uint16) for i in range(2)]
        if self.clock_type == 'slow clock':
            times = self.parent_device.change_times
        else:
            times = self.parent_device.times
        out_table = zeros(len(times),dtype=dtypes)
        out_table['freq0'].fill(1)
        out_table['freq1'].fill(1)
        
        for connection, dds in DDSs.items():
            out_table['freq%d'%connection] = dds.frequency.raw_output
            out_table['amp%d'%connection] = dds.amplitude.raw_output
            out_table['phase%d'%connection] = dds.phase.raw_output
        grp = hdf5_file.create_group('/devices/'+self.name)
        grp.attrs['frequency_scale_factor'] = 1e7
        grp.attrs['amplitude_scale_factor'] = 1023
        grp.attrs['phase_scale_factor'] = 45.511111111111113
        grp.create_dataset('TABLE_DATA',compression=compression,data=out_table) 

def generate_connection_table():
    all_devices = []
    connection_table = []
    devicedict = {}
    def sortkey(row):
        device = devicedict[row[0]]
        return str(device.generation) + device.name
    
    for device in inventory:
        devicedict[device.name] = device
        all_devices.extend(device.get_all_children())
        # The three child 'devices' of a DDS aren't really devices, ignore them.
        if not isinstance(device.parent_device,DDS):
            connection_table.append((device.name, device.__class__.__name__,
                                     device.parent_device.name if device.parent_device else str(None),
                                     str(device.connection if device.parent_device else str(None))))
    connection_table.sort(key=sortkey)
    connection_table_dtypes = [('name','a256'), ('class','a256'), ('parent','a256'), ('connected to','a256')]
    connection_table_array = empty(len(connection_table),dtype=connection_table_dtypes)
    for i, row in enumerate(connection_table):
        connection_table_array[i] = row
    hdf5_file.create_dataset('connection table',compression=compression,data=connection_table_array)
    print 'Name'.rjust(15), 'Class'.rjust(15), 'Parent'.rjust(15), 'connected to'.rjust(15)
    print '----'.rjust(15), '-----'.rjust(15), '------'.rjust(15), '------------'.rjust(15)
    for row in connection_table:
        print row[0].rjust(15), row[1].rjust(15), row[2].rjust(15), row[3].rjust(15)
                
def generate_code():
    for device in inventory:
        if not device.parent_device:
            device.generate_code()
            print
            print device.name + ':'
            print 'Fast clock'.ljust(15) + str(len(device.times)).rjust(15) + ' x ', str(device.times.dtype).ljust(15)
            print 'Slow clock'.ljust(15) + str(len(device.change_times)).rjust(15) + ' x ', str(device.change_times.dtype).ljust(15)
            for output in device.get_all_outputs():
                print output.name.ljust(15) + str(len(output.raw_output)).rjust(15) + ' x ', str(output.raw_output.dtype).ljust(15)
            print

def stop(t):
    for device in inventory:
        if not device.parent_device:  
            device.stop_time = t
    hdf5_file.create_group('/devices')
    generate_code()
    generate_connection_table()
    labscriptfile = os.path.join(sys.path[0],sys.argv[0])
    script = hdf5_file.create_dataset('script',compression=compression,data=open(labscriptfile).read())
    script.attrs['name'] = os.path.basename(sys.argv[0])
    script.attrs['path'] = sys.path[0]
    hdf5_file.close()
    
def open_hdf5_file():
    try:
        assert len(sys.argv) > 1
        hdf5_filename = sys.argv[-1]
        assert hdf5_filename.lower().endswith('h5')
    except:
        if not sys.path[0]:
            sys.stderr.write('ERROR: Can\'t run labscript interactively. Labscript relies on there being a script file. If you\'re just checking that you can import labscript, yes you can :)\n')
            sys.exit(1)
        newh5file = os.path.join(sys.path[0],sys.argv[0].split('.py')[0]+'.h5')
        if os.path.exists(newh5file) and '-replace' not in sys.argv:
            dialog = gtk.MessageDialog(None,0,gtk.MESSAGE_WARNING, gtk.BUTTONS_OK_CANCEL,
             'Replace existing hdf5 file ' + sys.argv[0].split('.py')[0]+'.h5?')
            dialog.present()
            response = dialog.run()
            # There is no gtk mainloop to make the dialog box
            # vanish. We'll do a little mainloop here until there are
            # no more events pending:
            dialog.destroy()
            while gtk.events_pending():
                gtk.main_iteration()
            if not response == gtk.RESPONSE_OK:
                sys.stderr.write('ERROR: No hdf5 file provided, and not overwriting previously existing h5 file with default filename. Stopping.\n')
                sys.exit(1)
        f = h5py.File(newh5file,'w')
        group = f.create_group('globals')
        f.close()
        hdf5_filename = newh5file
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
params = dict(hdf5_file['globals'].attrs)
if '-compress' in sys.argv:
    compression = 'gzip'
else:
    compression = None
        
for name in params.keys():
    if name in globals().keys() or name in locals().keys() or name in dir(__builtins__):
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         ' is already a name used by Python, labscript, or Pylab.'+
                         ' Please choose a different variable name to avoid a conflict.\n')
        sys.exit(1)
    if name in keyword.kwlist:
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         ' is a reserved Python keyword.' +
                         ' Please choose a different variable name.\n')
        sys.exit(1)
        
    try:
        assert '.' not in name
        exec(name + ' = 0')
        exec('del ' + name )
    except:
        sys.stderr.write('ERROR whilst parsing globals from %s. \'%s\''%(sys.argv[1],name) +
                         'is not a valid Python variable name.' +
                         ' Please choose a different variable name.\n')
        sys.exit(1)
    setattr(__main__,name, params[name])
    
if params.keys():
    # get rid of the loop variable -- it caused a subtle bug once by
    # continuing to exist:
    del name

       
