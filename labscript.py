import os
import sys
import subprocess
import keyword

import h5py
from pylab import *

import functions

import __main__

ns = 1e-9
us = 1e-6
ms = 1e-3
s = 1
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

# We need to backup the builtins as they are now, as well as have a
# reference to the actual builtins dictionary (which will change as we
# add globals and devices to it), so that we can restore the builtins
# when labscript_cleanup() is called. The if statement here is because
# __builtins__ is either a dictionary, or a module, depending on whether
# we are the __main__ module or not. We're usually not, but just in case
# we'll handle both cases.
module = type(__main__)
if isinstance(__builtins__, module):
    _builtins_dict = __builtins__.__dict__
    _existing_builtins_dict = __builtins__.__dict__.copy()
elif isinstance(__builtins__, dict):
    _builtins_dict = __builtins__
    _existing_builtins_dict = __builtins__.copy()
else:
    raise ImportError('Can\'t get __builtins__ dictionary, what\'s going on?')
    
# Startupinfo, for ensuring subprocesses don't launch with a visible command window:
if os.name=='nt':
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= 1 #subprocess.STARTF_USESHOWWINDOW # This variable isn't defined, but apparently it's equal to one.
else:
    startupinfo = None
        
        
class config:
    supress_mild_warnings = True
    supress_all_warnings = False
    compression = None # set to 'gzip' for compression 
   
    
class NoWarnings(object):
    """A context manager which sets config.supress_mild_warnings to True
    whilst in use.  Allows the user to supress warnings for specific
    lines when they know that the warning does not indicate a problem."""
    def __enter__(self):
        self.existing_warning_setting = config.supress_all_warnings
        config.supress_all_warnings = True
    def __exit__(self, *args):
        config.supress_all_warnings = self.existing_warning_setting
    
no_warnings = NoWarnings() # This is the object that should be used, not the class above


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
        compiler.inventory.append(self)
        if parent_device:
            parent_device.add_device(self)
        
        # Check that the name doesn't already exist in the python namespace
        if name in locals() or name in globals() or name in __builtins__:
            raise LabscriptError('The device name %s already exists in the Python namespace. Please choose another.'%name)
        if name in keyword.kwlist:
            raise LabscriptError('%s is a reserved Python keyword.'%name +
                                 ' Please choose a different device name.')
                                     
        try:
            # Test that name is a valid Python variable name:
            exec '%s = None'%name
            assert '.' not in name
        except:
            raise ValueError('%s is not a valid Python variable name.'%name)
        
        # Put self into the global namespace:
        __builtins__[name] = self
        
    def add_device(self,device):
        if any([isinstance(device,DeviceClass) for DeviceClass in self.allowed_children]):
            self.child_devices.append(device)
        else:
            raise LabscriptError('Devices of type %s cannot be attached to devices of type %s.'%(device.description,self.description))
            
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

    def generate_code(self, hdf5_file):
        for device in self.child_devices:
            device.generate_code(hdf5_file)
  

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
                raise LabscriptError('Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(change_times[i+1])) +
                                     'One or more connected devices cannot support update delays shorter than %s sec.'%str(1.0/self.clock_limit))
        # If the device has no children, we still need it to have a
        # single instruction. So we'll add 0 as a change time:
        if not change_times:
            change_times.append(0)
        # Also add the stop time as as change time. First check that it isn't too close to the time of the last instruction:
        dt = self.stop_time - change_times[-1]
        if abs(dt) < 1.0/self.clock_limit:
            raise LabscriptError('The stop time of the experiment is t= %s s, but the last instruction for a device attached to %s is at t= %s s. '%( str(self.stop_time), self.name, str(change_times[-1])) +
                                 'One or more connected devices cannot support update delays shorter than %s sec. Please set the stop_time a bit later.'%str(1.0/self.clock_limit))
        change_times.append(self.stop_time)
        # Sort change times so self.stop_time will be in the middle
        # somewhere if it is prior to the last actual instruction. Whilst
        # this means the user has set stop_time in error, not catching
        # the error here allows it to be caught later by the specific
        # device that has more instructions after self.stop_time. Thus
        # we provide the user with sligtly more detailed error info.
        change_times.sort()
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
                if hasattr(output,'timeseries') and isinstance(output.timeseries[i],dict) and output.timeseries[i]['clock rate'] > maxrate:
                    # It does have the highest clock rate? Then store that rate to max_rate:
                    maxrate = output.timeseries[i]['clock rate']
            if maxrate > self.clock_limit:
                raise LabscriptError('At t = %s sec, a clock rate of %s Hz was requested. '%(str(time),str(maxrate)) + 
                                    'One or more devices connected to %s cannot support clock rates higher than %sHz.'%(str(self.name),str(self.clock_limit)))
                
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
                all_times.append(ticks)
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
                    if self.stop_time > time:
                        # There is no next instruction. Hold the last clock
                        # tick until self.stop_time.
                        raise Exception('This shouldn\'t happen -- stop_time should always be equal to the time of the last instruction. Please report a bug.')
                        clock.append({'start': time, 'reps': 1, 'step': self.stop_time - time,'slow_clock_tick':True})
                    # Error if self.stop_time has been set to less
                    # than the time of the last instruction:
                    elif self.stop_time < time:
                        raise LabscriptError('%s %s has more instructions after the experiment\'s stop time.'%(self.description,self.name))
                    # If self.stop_time is the same as the time of the last
                    # instruction, then we'll get the last instruction
                    # out still, so that the total number of clock
                    # ticks matches the number of data points in the
                    # Output.raw_output arrays. We'll make this last
                    # cycle be at ten times the minimum step duration.
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
        self.change_times = fastflatten(change_times, float)
        self.times = fastflatten(all_times,float)
        
    def generate_code(self, hdf5_file):
        self.generate_clock()
        Device.generate_code(self, hdf5_file)
        
        
class PulseBlaster(PseudoClock):
    
    pb_instructions = {'CONTINUE':   0,
                       'STOP':       1, 
                       'LOOP':       2, 
                       'END_LOOP':   3,
                       'BRANCH':     6,
                       'LONG_DELAY': 7,
                       'WAIT':       8}
                       
    description = 'PB-DDSII-300'
    clock_limit = 8.3e6 # Slight underestimate I think.
    fast_clock_flag = 0
    slow_clock_flag = 1
    clock_type = 'slow clock'
    
    # This value is coupled to a value in the PulseBlaster worker process of BLACS
    # This number was found experimentally but is determined theoretically by the
    # instruction lengths in BLACS, and a finite delay in the PulseBlaster
    #
    # IF YOU CHANGE ONE, YOU MUST CHANGE THE OTHER!
    trigger_delay = 250e-9 
    
    def __init__(self,name,board_number):
        PseudoClock.__init__(self,name,None,None)
        self.BLACS_connection = board_number
        
        self.trigger_data = {}
    
    def trigger(self,t,channel):
        if self.trigger_data:
            raise LabscriptError('You can only trigger the PulseBlaster %s once. If you require more than one trigger, please file a feature request for experiment WAITS on redmine'%self.name)
    
        if channel == 'software':    
            if t != 0:
                raise LabscriptError('You can only software trigger %s at t=0'%self.name)
            else:
                self.trigger_data[t] = {'type':'software'}
            
            return 0 # software triggers take no time
        else:
            # if we are not triggering at the start, make sure the channel stays high otherwise it will be set to 0 as default, which would be bad
            if t != 0:
                channel.go_high(0)
                # Just make sure that the time to trigger is not
                if t < 1.0/channel.clock_limit:
                    raise LabscriptError('The trigger %s used to trigger %s is too close to t=0 to change state. Please trigger %s at either t=0 or t>=%s'%(channel.name,self.name,self.name,1.0/channel.clock_limit))
            channel.go_low(t)
            channel.go_high(t+1.0/channel.clock_limit)
            self.trigger_data[t] = {'type':'hardware', 'channel':channel.name, 'device':channel.parent_device.name}
            
            return 2.0/channel.clock_limit # return the time this trigger will take (2 clock ticks)
    
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
                    raise LabscriptError('%s %s has invalid connection string: \'%s\'. '%(output.description,output.name,str(output.connection)) + 
                                         'Format must be \'flag n\' with n an integer less than 12, or \'dds n\' with n less than 2.')
                if not connection < 12:
                    raise LabscriptError('%s is set as connected to output connection %d of %s. '%(output.name, connection, self.name) +
                                         'Output connection number must be a integer less than 12.')
                if prefix == 'dds' and not connection < 2:
                    raise LabscriptError('%s is set as connected to output connection %d of %s. '%(output.name, connection, self.name) +
                                         'DDS output connection number must be a integer less than 2.')
                if prefix == 'flag' and connection in [self.slow_clock_flag, self.fast_clock_flag]:
                    raise LabscriptError('%s is set as connected to flag %d of %s.'%(output.name, connection, self.name) +
                                         'This is one of the PulseBlaster\'s clock flags.')
                for other_output in dig_outputs + dds_outputs:
                    if output.connection == other_output.connection:
                        raise LabscriptError('%s and %s are both set as connected to %s of %s.'%(output.name, other_output.name, output.connection, self.name))
                if isinstance(output,DigitalOut):
                	dig_outputs.append(output)
                elif isinstance(output, DDS):
                	dds_outputs.append(output)
                
        return dig_outputs, dds_outputs

    def generate_registers(self, hdf5_file, dds_outputs):
        ampdicts = {}
        phasedicts = {}
        freqdicts = {}
        group = hdf5_file.create_group('/devices/'+self.name)
        dds_dict = {}
        for output in dds_outputs:
            num = int(output.connection.split()[1])
            dds_dict[num] = output
        for num in [0,1]:
            
            if num in dds_dict:
                output = dds_dict[num]
            
                # Ensure that amplitudes are within bounds:
                if any(output.amplitude.raw_output > 1)  or any(output.amplitude.raw_output < 0):
                    raise LabscriptError('%s %s '%(output.amplitude.description, output.amplitude.name) +
                                      'can only have values between 0 and 1, ' + 
                                      'the limit imposed by %s.'%output.name)
                                      
                # Ensure that frequencies are within bounds:
                if any(output.frequency.raw_output > 150e6 )  or any(output.frequency.raw_output < 0):
                    raise LabscriptError('%s %s '%(output.frequency.description, output.frequency.name) +
                                      'can only have values between 0Hz and and 150MHz, ' + 
                                      'the limit imposed by %s.'%output.name)
                                      
                # Ensure that phase wraps around:
                output.phase.raw_output %= 360
                
                amps = set(output.amplitude.raw_output)
                phases = set(output.phase.raw_output)
                freqs = set(output.frequency.raw_output)
            else:
                # If the DDS is unused, it will use the following values
                # for the whole experimental run:
                amps = set([0])
                phases = set([0])
                freqs = set([0])
                                  
            if len(amps) > 1024:
                raise LabscriptError('%s dds%d can only support 1024 amplitude registers, and %s have been requested.'%(self.name, num, str(len(amps))))
            if len(phases) > 128:
                raise LabscriptError('%s dds%d can only support 128 phase registers, and %s have been requested.'%(self.name, num, str(len(phases))))
            if len(freqs) > 1024:
                raise LabscriptError('%s dds%d can only support 1024 frequency registers, and %s have been requested.'%(self.name, num, str(len(freqs))))
                                
            # start counting at 1 to leave room for the dummy instruction,
            # which BLACS will fill in with the state of the front
            # panel:
            ampregs = range(1,len(amps)+1)
            freqregs = range(1,len(freqs)+1)
            phaseregs = range(1,len(phases)+1)
            
            ampdicts[num] = dict(zip(amps,ampregs))
            freqdicts[num] = dict(zip(freqs,freqregs))
            phasedicts[num] = dict(zip(phases,phaseregs))
            
            # The zeros are the dummy instructions:
            freq_table = array([0] + list(freqs), dtype = float64) / 1e6 # convert to MHz
            amp_table = array([0] + list(amps), dtype = float32)
            phase_table = array([0] + list(phases), dtype = float64)
            
            subgroup = group.create_group('DDS%d'%num)
            subgroup.create_dataset('FREQ_REGS', compression=config.compression,data = freq_table)
            subgroup.create_dataset('AMP_REGS', compression=config.compression,data = amp_table)
            subgroup.create_dataset('PHASE_REGS', compression=config.compression,data = phase_table)
            
        return freqdicts, ampdicts, phasedicts
        
    def convert_to_pb_inst(self, hdf5_file, dig_outputs, dds_outputs, freqs, amps, phases):
        pb_inst = []
        # An array for storing the line numbers of the instructions at
        # which the slow clock ticks:
        slow_clock_indices = []
        # index to keep track of where in output.raw_output the
        # pulseblaster flags are coming from
        i = 0
        # index to record what line number of the pulseblaster hardware
        # instructions we're up to:
        j = 0
        # We've delegated the initial two instructions off to BLACS, which
        # can ensure continuity with the state of the front panel. Thus
        # these two instructions don't actually do anything:
        flags = [0]*12
        freqregs = [0]*2
        ampregs = [0]*2
        phaseregs = [0]*2
        dds_enables = [0]*2
        flags[self.fast_clock_flag] = 0
        flags[self.slow_clock_flag] = 0 
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
                        'flags': ''.join([str(flag) for flag in flags]), 'instruction': 'STOP',
                        'data': 0, 'delay': 10.0/self.clock_limit*1e9})
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
                        'flags': ''.join([str(flag) for flag in flags]), 'instruction': 'STOP',
                        'data': 0, 'delay': 10.0/self.clock_limit*1e9})    
        j += 2
        flagstring = '000000000000' # So that this variable is still defined if the for loop has no iterations
        for k, instruction in enumerate(self.clock):
            flags = [0]*12
            # The registers below are ones, not zeros, so that we don't
            # use the BLACS-inserted initial instructions. Instead
            # unused DDSs have a 'zero' in register one for freq, amp
            # and phase.
            freqregs = [1]*2
            ampregs = [1]*2
            phaseregs = [1]*2
            dds_enables = [0]*2
            for output in dig_outputs:
                flagindex = int(output.connection.split()[1])
                flags[flagindex] = int(output.raw_output[i])
            for output in dds_outputs:
                ddsnumber = int(output.connection.split()[1])
                freqregs[ddsnumber] = freqs[ddsnumber][output.frequency.raw_output[i]]
                ampregs[ddsnumber] = amps[ddsnumber][output.amplitude.raw_output[i]]
                phaseregs[ddsnumber] = phases[ddsnumber][output.phase.raw_output[i]]
                dds_enables[ddsnumber] = output.gate.raw_output[i]
            flags[self.fast_clock_flag] = 1
            flags[self.slow_clock_flag] = 1 if instruction['slow_clock_tick'] else 0
            if instruction['slow_clock_tick']:
                slow_clock_indices.append(j)
            flagstring = ''.join([str(flag) for flag in flags])
            if instruction['reps'] > 1048576:
                raise LabscriptError('Pulseblaster cannot support more than 1048576 loop iterations. ' +
                                      str(instruction['reps']) +' were requested at t = ' + str(instruction['start']) + '. '+
                                     'This can be fixed easily enough by using nested loops. If it is needed, ' +
                                     'please file a feature request at' +
                                     'http://redmine.physics.monash.edu.au/projects/labscript.')
                
            # Instruction delays > 55 secs will require a LONG_DELAY
            # to be inserted. How many times does the delay of the
            # loop/endloop instructions go into 55 secs?
            quotient, remainder = divmod(instruction['step']/2.0,55.0)
            if remainder < 100e-9:
                # The remainder will be used for the total duration of the LOOP and END_LOOP instructions. 
                # It must not be too short for this, if it is, take one LONG_DELAY iteration and give 
                # its duration to the loop instructions:
                quotient, remainder = quotient - 1, remainder + 55.0
                assert quotient >= 0 # Something is wrong if this is not the case
                
            # The loop and endloop instructions will only use the remainder:
            pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
                            'flags': flagstring, 'instruction': 'LOOP',
                            'data': instruction['reps'], 'delay': remainder*1e9})
            flags[self.fast_clock_flag] = 0
            flags[self.slow_clock_flag] = 0
            flagstring = ''.join([str(flag) for flag in flags])
            # If there was a nonzero quotient, let's wait twice that
            # many multiples of 55 seconds (one multiple of 55 seconds
            # for each of the other two loop and endloop instructions):
            if quotient:
                pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
                            'flags': flagstring, 'instruction': 'LONG_DELAY',
                            'data': int(2*quotient), 'delay': 55*1e9})
                            
            pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
                            'flags': flagstring, 'instruction': 'END_LOOP',
                            'data': j, 'delay': remainder*1e9})
                            
            # Two instructions were used in the case of there being no LONG_DELAY, 
            # otherwise three. This increment is done here so that the j referred
            # to in the previous line still refers to the LOOP instruction.
            j += 3 if quotient else 2
            try:
                if self.clock[k+1]['slow_clock_tick']:
                    i += 1
            except IndexError:
                pass
        # This is how we stop the pulse program. We branch from the last
        # instruction to the zeroth, which BLACS has programmed in with
        # the same values and a WAIT instruction. The PulseBlaster then
        # waits on instuction zero, which is a state ready for either
        # further static updates or buffered mode.
        pb_inst.append({'freqs': freqregs, 'amps': ampregs, 'phases': phaseregs, 'enables':dds_enables,
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
            en0 = inst['enables'][0]
            en1 = inst['enables'][1]
            pb_inst_table[i] = (freq0,phase0,amp0,en0,0,freq1,phase1,amp1,en1,0, flagint, 
                                instructionint, dataint, delaydouble)
        slow_clock_indices = array(slow_clock_indices, dtype = uint32)                  
        # Okey now write it to the file: 
        group = hdf5_file['/devices/'+self.name]  
        group.create_dataset('PULSE_PROGRAM', compression=config.compression,data = pb_inst_table)         
        group.create_dataset('FAST_CLOCK', compression=config.compression,data = self.times)         
        group.create_dataset('SLOW_CLOCK', compression=config.compression,data = self.change_times)   
        group.create_dataset('CLOCK_INDICES', compression=config.compression,data = slow_clock_indices)  
        group.attrs['stop_time'] = self.stop_time       
                              
    def generate_code(self, hdf5_file):
        # since we only have one trigger point at the moment (because WAITS are not implemented)
        # let's get the trigger time:
        if not self.trigger_data:
            raise LabscriptError('%s %s has not been triggered. Please provide a trigger instruction.'%(self.description, self.name))
        first_trigger_time, first_trigger_data = self.trigger_data.items()[0]
        
        #eventually this will be put in a loop, that iterates over the entries in the trigger_data dictionary!
        # We add on the trigger delay here to remove the time offset introduced by the length of the WAIT
        # instructions BLACS adds at the start of the pulse program, and the finite trigger time inherent in the pulseBlaster
        trigger_time = first_trigger_time + (self.trigger_delay if first_trigger_data['type'] != 'software' else 0)
        # This is stored so that child devices know what time to add an initial instruction in at, if none was specified by the user
        self.trigger_time = trigger_time
    
        # Get all outputs (on pseudo clock and direct)        
        # Iterate over outputs
        for output in self.get_all_outputs():
            offset_instructions = {}
            for time,instruction in output.instructions.items():            
                # Check that no instruction is before the first trigger
                if time < (first_trigger_time + (self.trigger_delay if first_trigger_data['type'] != 'software' else 0)):
                    if time < first_trigger_time:   
                        raise LabscriptError('You cannot command output %s (at t=%s) before triggering it\'s clock (device %s at t=%s)'%(output.name,str(time),self.name,str(first_trigger_time+self.trigger_delay)))
                    else:
                        raise LabscriptError('Triggering %s takes %ss. Thus, you cannot command %s before t=%ss (you have attempted to command it at t=%ss)'%(self.name,str(self.trigger_delay),output.name,str(first_trigger_time+self.trigger_delay),str(time)))
                
                        
                # check that no loop bridges a trigger
                if isinstance(instruction,dict):
                    if time < trigger_time and instruction['end time'] >= trigger_time:
                        raise LabscriptError('You cannot trigger %s at t=%s during a ramp from t=%s to t=%s on %s'%(self.name,str(trigger_time),str(time),str(instruction['end time']),output.name))
                
                # check that no instruction is too close to the trigger (both before and after)
                
                
                # Offset the times in the instruction list from the trigger.
                if isinstance(instruction,dict):
                    # Note that this will actually modify output.instructions. Not that it matters,
                    # But if you change the code, it might!
                    instruction['end time'] = instruction['end time'] - trigger_time
                    instruction['initial time'] = instruction['initial time'] - trigger_time
                
                offset_instructions[time-trigger_time] = instruction
            output.instructions = offset_instructions
    
        # Adjust the stop time relative to the last trigger time
        self.stop_time = self.stop_time-trigger_time
    
        # Generate the hardware instructions
        PseudoClock.generate_code(self, hdf5_file)
        dig_outputs, dds_outputs = self.get_direct_outputs()
        freqs, amps, phases = self.generate_registers(hdf5_file, dds_outputs)
        self.convert_to_pb_inst(hdf5_file, dig_outputs, dds_outputs, freqs, amps, phases)
        
        # Make sure the device has been triggered, and save the data in the HDF5 file
        if self.trigger_data:
            group = hdf5_file['/devices/'+self.name]
            group.attrs['trigger'] = repr(self.trigger_data)
        else:
            raise LabscriptError('You must trigger %s sometime during the experiment'%self.name)
        

            
class Output(Device):
    description = 'generic output'
    allowed_states = {}
    dtype = float32
    scale_factor = 1
    generation = 3
    def __init__(self,name,parent_device,connection,limits = None,unit_conversion_class = None,unit_conversion_parameters = None):
        self.instructions = {}
        self.ramp_limits = [] # For checking ramps don't overlap
        self.clock_type = parent_device.clock_type
        if not unit_conversion_parameters:
            unit_conversion_parameters = {}
        self.unit_conversion_class = unit_conversion_class
        self.unit_conversion_parameters = unit_conversion_parameters        
        Device.__init__(self,name,parent_device,connection)  
        
        # Instatiate the calibration
        if unit_conversion_class is not None:
            self.calibration = unit_conversion_class(unit_conversion_parameters)
            # Validate the calibration class
            for units in self.calibration.derived_units:
                #Does the conversion to base units function exist for each defined unit type?
                if not hasattr(self.calibration,units+"_to_base"):
                    raise LabscriptError('The function "%s_to_base" does not exist within the calibration "%s" used in output "%s"'%(units,self.unit_conversion_class,self.name))
                #Does the conversion to base units function exist for each defined unit type?
                if not hasattr(self.calibration,units+"_from_base"):
                    raise LabscriptError('The function "%s_from_base" does not exist within the calibration "%s" used in output "%s"'%(units,self.unit_conversion_class,self.name))
        
        # If limits exist, check they are valid
        # Here we specifically differentiate "None" from False as we will later have a conditional which relies on
        # self.limits being either a correct tuple, or "None"
        if limits is not None:
            if not isinstance(limits,tuple) or len(limits) is not 2:
                raise LabscriptError('The limits for "%s" must be tuple of length 2. Eg. limits=(1,2)'%(self.name))
            if limits[0] > limits[1]:
                raise LabscriptError('The first element of the tuple must be lower than the second element. Eg limits=(1,2), NOT limits=(2,1)')
        # Save limits even if they are None        
        self.limits = limits
    
    def triger(self,t,*args,**kwargs):
        raise LabscriptError('You cannot trigger the channel %s at time %s. Did you mean to trigger its parent %s?'%(self.name,str(t),self.parent_device.name))
    
    @property
    def clock_limit(self):
        parent = self.find_pseudoclock()
        return parent.clock_limit
    
    def find_pseudoclock(self):
        parent = self.parent_device
        while parent.parent_device != None:
            parent = parent.parent_device
            
        if isinstance(parent,PseudoClock):
            return parent
        else:
            raise LabscriptError('The parent of %s is not a pseudoclock. Labscript requires that toplevel devices (with no parent) are a pseudoclock.'%self.name)
    
    
    def apply_calibration(self,value,units):
        # Is a calibration in use?
        if self.unit_conversion_class is None:
            raise LabscriptError('You can not specify the units in an instruction for output "%s" as it does not have a calibration associated with it'%(self.name))
                    
        # Does a calibration exist for the units specified?
        if units not in self.calibration.derived_units:
            raise LabscriptError('The units "%s" does not exist within the calibration "%s" used in output "%s"'%(units,self.unit_conversion_class,self.name))
                    
        # Return the calibrated value
        return getattr(self.calibration,units+"_to_base")(value)
        
    def instruction_to_string(self,instruction):
        """gets a human readable description of an instruction"""
        if isinstance(instruction,dict):
            return instruction['description']
        elif self.allowed_states:
            return str(self.allowed_states[instruction])
        else:
            return str(instruction)

    def add_instruction(self,time,instruction,units=None):
        # round to the nearest 0.1 nanoseconds, to prevent floating point
        # rounding errors from breaking our equality checks later on.
        time = round(time,10)
        # Also round end time of ramps to the nearest 0.1 ns:
        if isinstance(instruction,dict):
            instruction['end time'] = round(instruction['end time'],10)
            instruction['initial time'] = round(instruction['initial time'],10)
        #Check that time is not negative:
        if time < 0:
            err = ' '.join([self.description, self.name, 'has an instruction at t=%ss.'%str(time),
                 'Labscript cannot handle instructions earlier than the experiment start time (t=0).'])
            raise LabscriptError(err)
        #Check that this doesn't collide with previous instructions:
        if time in self.instructions.keys():
            if not config.supress_all_warnings:
                message = ' '.join(['WARNING: State of', self.description, self.name, 'at t=%ss'%str(time),
                          'has already been set to %s.'%self.instruction_to_string(self.instructions[time]),
                          'Overwriting to %s. (note: all values in base units where relevant)'%self.instruction_to_string(self.apply_calibration(instruction,units) if units and not isinstance(instruction,dict) else instruction)])
                sys.stderr.write(message+'\n')
        # Check that ramps don't collide
        if isinstance(instruction,dict):
            # No ramps allowed if this output is on a slow clock:
            if self.clock_type == 'slow clock':
                raise LabscriptError('%s %s is on a slow clock.'%(self.description, self.name) + 
                                     'It cannot have a function ramp as an instruction.')
            for start, end in self.ramp_limits:
                if start < time < end or start < instruction['end time'] < end:
                    err = ' '.join(['State of', self.description, self.name, 'from t = %ss to %ss'%(str(start),str(end)),
                        'has already been set to %s.'%self.instruction_to_string(self.instructions[start]),
                        'Cannot set to %s from t = %ss to %ss.'%(self.instruction_to_string(instruction),str(time),str(instruction['end time']))])
                    raise LabscriptError(err)
                self.ramp_limits.append((time,instruction['end time']))
            # Check that start time is before end time:
            if time > instruction['end time']:
                raise LabscriptError('%s %s has been passed a function ramp %s with a negative duration.'%(self.description, self.name, self.instruction_to_string(instruction)))
            if instruction['clock rate'] == 0:
                raise LabscriptError('A nonzero sample rate is required.')
            # Else we have a "constant", single valued instruction
        else:
            # If we have units specified, convert the value
            if units is not None:
                # Apply the unit calibration now
                instruction = self.apply_calibration(instruction,units)
            # if we have limits, check the value is valid
            if self.limits:
                if self.limits[0] <= instruction <= self.limits[1]:
                    raise LabscriptError('You cannot program the value %d (base units) to %s as it falls outside the limits (%d to %d)'%(instruction,self.name,limits[0],limits[1]))
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
            if not config.supress_mild_warnings and not config.supress_all_warnings:
                sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions. It will be set to %s for all time.\n'%self.instruction_to_string(self.default_value)]))
            self.add_instruction(0,self.default_value)  
        # Check if there are no instructions at t=0. Generate a warning and insert an
        # instruction telling the output to start at zero.
        if 0 not in self.instructions.keys():
            if not config.supress_mild_warnings and not config.supress_all_warnings:
               sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions at t=0. It will initially be set to %s.\n'%self.instruction_to_string(0)]))
            self.add_instruction(0,self.default_value) 
        # Check that ramps have instructions following them.
        # If they don't, insert an instruction telling them to hold their final value.
        for instruction in self.instructions.values():
            if isinstance(instruction, dict) and instruction['end time'] not in self.instructions.keys():
                self.add_instruction(instruction['end time'], instruction['function'](instruction['end time']-instruction['initial time']),instruction['units'])
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
                    outarray = self.timeseries[i]['function'](midpoints-self.timeseries[i]['initial time'])
                    # Now that we have the list of output points, pass them through the unit calibration
                    if self.timeseries[i]['units'] is not None:
                        outarray = self.apply_calibration(outarray,self.timeseries[i]['units'])
                    # if we have limits, check the value is valid
                    if self.limits:
                        if ((outarray<self.limits[0])|(outarray>self.limits[1])).any():
                            raise LabscriptError('The function %s called on "%s" at t=%d generated a value which falls outside the base unit limits (%d to %d)'%(self.timeseries[i]['function'],self.name,midpoints[0],limits[0],limits[1]))
                else:
                    outarray = empty(len(time),dtype=float32)
                    outarray.fill(self.timeseries[i])
                outputarray.append(outarray)
            else:
                outputarray.append(self.timeseries[i])
        del self.timeseries # don't need this any more.
        self.raw_output = fastflatten(outputarray, self.dtype)
        

class AnalogQuantity(Output):
    description = 'analog quantity'
    default_value = 0
    def ramp(self,t,duration,initial,final,samplerate,units=None):
        self.add_instruction(t, {'function': functions.ramp(duration,initial,final), 'description':'linear ramp',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})
        
        return duration
                                 
    def sine(self,t,duration,amplitude,angfreq,phase,dc_offset,samplerate,units=None):
        self.add_instruction(t, {'function': functions.sine(duration,amplitude,angfreq,phase,dc_offset), 'description':'sine wave',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})
       
        return duration
        
    def sine_ramp(self,t,duration,initial,final,samplerate,units=None):
        self.add_instruction(t, {'function': functions.sine_ramp(duration,initial,final), 'description':'sinusoidal ramp',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})   
                
        return duration
    
    def exp_ramp(self, t, duration, initial, final, samplerate, zero=0, trunc=False, trunc_type='linear', units=None):
        if trunc is not False:
            if trunc_type == 'linear':
                trunc_duration = duration*log((initial-zero)/(trunc-zero))/log((initial-zero)/(final-zero))
            if trunc_type == 'exponential':
                trunc_duration = trunc * duration
                # final = functions.exp_ramp(0, duration, initial, final, zero)(trunc_duration)
            else:
                raise LabscriptError('Truncation type for exp_ramp not supported. Must be either linear or exponential.')
        else:
            trunc_duration = duration
        self.add_instruction(t, {'function': functions.exp_ramp(duration,initial,final,zero), 'description':'exponential ramp',
                             'initial time':t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
        
        return trunc_duration
     
    def exp_ramp_t(self, t, duration, initial, final, time_constant, samplerate, trunc=False, units=None):
        # Exponential ramp set by the time constant. No truncation yet!
        zero = (final-initial*exp(-duration/time_constant)) / (1-exp(-duration/time_constant))
        if trunc is not False:
            trunc_duration = time_constant * log((initial-zero)/(trunc-zero))
        else:
            trunc_duration = duration
        self.add_instruction(t, {'function': functions.exp_ramp_t(duration, initial, final, time_constant), 'description':'exponential ramp with time consntant',
                             'initial time':t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
                
        return trunc_duration
                             
    def constant(self,t,value,units=None):
        # verify that value can be converted to float
        val = float(value)
        self.add_instruction(t,value,units)
        
      
class AnalogOut(AnalogQuantity):
    description = 'analog output'
    
class StaticAnalogQuantity(Output):
    description = 'static analog quantity'
    default_value = 0.0
    
    value_set = False
    
    def constant(self,value,units=None):
        if not self.value_set:
            # Since this is a StaticAnalogOut, we set the instruction to be at t=0.0
            #self.add_instruction(0.0,value,units)
            
            # If we have units specified, convert the value
            if units is not None:
                # Apply the unit calibration now
                value = self.apply_calibration(value,units)
            # if we have limits, check the value is valid
            if self.limits:
                if self.limits[0] <= value <= self.limits[1]:
                    raise LabscriptError('You cannot program the value %d (base units) to %s as it falls outside the limits (%d to %d)'%(value,self.name,limits[0],limits[1]))
        
            
            self.value_set = True
            self.static_value = value
        else:
            raise LabscriptError('%s %s has already been set to %s (base units). It cannot also be set to %s (%s).'%(self.description, self.name, str(self.static_value), str(value),units if units is not None else "base units"))
    
    def get_change_times(self):
        if not hasattr(self,'static_value'):
            self.static_value = self.default_value
            
        return [] # Return an empty list as the calling function at the pseudoclock level expects a list
    
    def make_timeseries(self,change_times):
        pass
    
    def expand_timeseries(self,*args,**kwargs):
        self.raw_output = array([self.static_value])

class StaticAnalogOut(StaticAnalogQuantity):
    description = 'static analog output'
        
class DigitalQuantity(Output):
    description = 'digital quantity'
    allowed_states = {1:'high', 0:'low'}
    default_value = 0
    dtype = uint32
    
    # Redefine __init__ so that you cannot define a limit of calibration for DO
    def __init__(self,name,parent_device,connection):
        Output.__init__(self,name,parent_device,connection)
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 

class DigitalOut(DigitalQuantity):
    description = 'digital output'

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
        return end_time - start_time
  
class IntermediateDevice(Device):
    generation = 1
    
    def __init__(self, name, parent_device,clock_type):
        self.name = name
        if not clock_type in ['fast clock', 'slow clock']:
            raise LabscriptError("Clock type for %s %s can only be 'slow clock' or 'fast clock'."%(self.name,self.description))
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
    
    def __init__(self, name, parent_device, clock_type, clock_terminal, acquisition_rate=0):
        IntermediateDevice.__init__(self, name, parent_device,clock_type)
        self.acquisition_rate = acquisition_rate
        self.clock_terminal = clock_terminal
        
    def add_device(self,output):
        # TODO: check there are no duplicates, check that connection
        # string is formatted correctly.
        Device.add_device(self,output)
        
    def convert_bools_to_bytes(self,digitals):
        """converts digital outputs to an array of bitfields stored
        as self.digital_dtype"""
        outputarray = [0]*self.n_digitals
        for output in digitals:
            port, line = output.connection.replace('port','').replace('line','').split('/')
            port, line  = int(port),int(line)
            if port > 0:
                raise LabscriptError('Ports > 0 on NI Boards not implemented. Please use port 0, or file a feature request at redmine.physics.monash.edu.au/labscript.')
            outputarray[line] = output.raw_output
        bits = bitfield(outputarray,dtype=self.digital_dtype)
        return bits
            
    def generate_code(self, hdf5_file):
        Device.generate_code(self, hdf5_file)
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
                raise LabscriptError('%s %s '%(output.description, output.name) +
                                  'can only have values between -10 and 10 Volts, ' + 
                                  'the limit imposed by %s.'%self.name)
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
        digital_out_table = []
        if digitals:
            digital_out_table = self.convert_bools_to_bytes(digitals.values())
        grp = hdf5_file.create_group('/devices/'+self.name)
        if all(analog_out_table.shape): # Both dimensions must be nonzero
            analog_dataset = grp.create_dataset('ANALOG_OUTS',compression=config.compression,data=analog_out_table)
            grp.attrs['analog_out_channels'] = ', '.join(analog_out_attrs)
        if len(digital_out_table): # Table must be non empty
            digital_dataset = grp.create_dataset('DIGITAL_OUTS',compression=config.compression,data=digital_out_table)
            grp.attrs['digital_lines'] = '/'.join((self.name,'port0','line0:%d'%(self.n_digitals-1)))
        if len(acquisition_table): # Table must be non empty
            input_dataset = grp.create_dataset('ACQUISITIONS',compression=config.compression,data=acquisition_table)
            grp.attrs['analog_in_channels'] = ', '.join(input_attrs)
            grp.attrs['acquisition_rate'] = self.acquisition_rate
        grp.attrs['clock_terminal'] = self.clock_terminal
        
        
class NI_PCI_6733(NIBoard):
    description = 'NI-PCI-6733'
    n_analogs = 8
    n_digitals = 0
    n_analog_ins = 0
    digital_dtype = uint32
    
    def generate_code(self, hdf5_file):
        NIBoard.generate_code(self, hdf5_file)
        if len(self.child_devices) % 2:
            raise LabscriptError('%s %s must have an even numer of analog outputs '%(self.description, self.name) +
                             'in order to guarantee an even total number of samples, which is a limitation of the DAQmx library. ' +
                             'Please add a dummy output device or remove an output you\'re not using, so that there are an even number of outputs. Sorry, this is annoying I know :).')
        
                        
class NI_PCIe_6363(NIBoard):
    description = 'NI-PCIe-6363'
    n_analogs = 4
    n_digitals = 32
    n_analog_ins = 32
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

    def generate_code(self, hdf5_file):
        classname = self.__class__.__name__
        calibration_table_dtypes = [('name','a256'), ('open_delay',float), ('close_delay',float)]
        if classname not in hdf5_file['calibrations']:
            hdf5_file['calibrations'].create_dataset(classname, (0,), dtype=calibration_table_dtypes, maxshape=(None,))
        metadata = (self.name,self.open_delay,self.close_delay)
        dataset = hdf5_file['calibrations'][classname]
        dataset.resize((len(dataset)+1,))
        dataset[len(dataset)-1] = metadata
        
        
class Camera(DigitalOut):
    description = 'Generic Camera'
    frame_types = ['atoms','flat','dark','fluoro','clean']
    minimum_recovery_time = 0 # To be set by subclasses
    
    def __init__(self,name,parent_device,connection,BIAS_port,serial_number,SDK,effective_pixel_size,exposuretime,orientation):
        DigitalOut.__init__(self,name,parent_device,connection)
        self.exposuretime = exposuretime
        self.orientation = orientation
        self.exposures = []
        self.go_low(0)
        self.BLACS_connection = BIAS_port
        if isinstance(serial_number,str):
            serial_number = int(serial_number,16)
        self.sn = uint64(serial_number)
        self.sdk = str(SDK)
        self.effective_pixel_size = effective_pixel_size
        
    def expose(self,name, t ,frametype):
        self.go_high(t)
        duration = self.exposuretime
        self.go_low(t + duration)
        for exposure in self.exposures:
            start = exposure[1]
            end = start + duration
            # Check for overlapping exposures:
            if start <= t <= end or start <= t+duration <= end:
                raise LabscriptError('%s %s has two overlapping exposures: ' %(self.description, self.name) + \
                                 'one at t = %fs for %fs, and another at t = %fs for %fs.'%(t,duration,start,duration))
            # Check for exposures too close together:
            if abs(start - (t + duration)) < self.minimum_recovery_time or abs((t+duration) - end) < self.minimum_recovery_time:
                raise LabscriptError('%s %s has two exposures closer together than the minimum recovery time: ' %(self.description, self.name) + \
                                 'one at t = %fs for %fs, and another at t = %fs for %fs. '%(t,duration,start,duration) + \
                                 'The minimum recovery time is %fs.'%self.minimum_recovery_time)
        # Check for invalid frame type:                        
        if not frametype in self.frame_types:
            raise LabscriptError('%s is not a valid frame type for %s %s.'%(str(frametype), self.description, self.name) +\
                             'Allowed frame types are: \n%s'%'\n'.join(self.frame_types))
        self.exposures.append((name, t, frametype))
        
    def generate_code(self, hdf5_file):
        table_dtypes = [('name','a256'), ('time',float), ('frametype','a256')]
        data = array(self.exposures,dtype=table_dtypes)
        group = hdf5_file['devices'].create_group(self.name)
        group.attrs['exposure_time'] = float(self.exposuretime)
        group.attrs['orientation'] = self.orientation
        group.attrs['SDK'] = self.sdk
        group.attrs['serial_number'] = self.sn
        group.attrs['effective_pixel_size'] = self.effective_pixel_size
        if self.exposures:
            group.create_dataset('EXPOSURES', data=data)

            
class DDS(Device):
    description = 'DDS'
    allowed_children = [AnalogQuantity,DigitalOut,DigitalQuantity] # Adds its own children when initialised
    generation = 2
    def __init__(self,name,parent_device,connection,digital_gate={},freq_limits = None,freq_conv_class = None,freq_conv_params = {},amp_limits=None,amp_conv_class = None,amp_conv_params = {},phase_limits=None,phase_conv_class = None,phase_conv_params = {}):
        self.clock_type = parent_device.clock_type
        Device.__init__(self,name,parent_device,connection)
        self.frequency = AnalogQuantity(self.name+'_freq',self,'freq',freq_limits,freq_conv_class,freq_conv_params)
        self.amplitude = AnalogQuantity(self.name+'_amp',self,'amp',amp_limits,amp_conv_class,amp_conv_params)
        self.phase = AnalogQuantity(self.name+'_phase',self,'phase',phase_limits,phase_conv_class,phase_conv_params)
        self.gate = None
        if isinstance(self.parent_device,NovaTechDDS9M):
            self.frequency.default_value = 0.1
            if 'device' in digital_gate and 'connection' in digital_gate:            
                self.gate = DigitalOut(self.name+'_gate',digital_gate['device'],digital_gate['connection'])
            # Did they only put one key in the dictionary, or use the wrong keywords?
            elif len(digital_gate) > 0:
                raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.'%(self.name))
        elif isinstance(self.parent_device,PulseBlaster):
            if 'device' in digital_gate and 'connection' in digital_gate: 
                raise LabscriptError('You cannot specify a digital gate for a DDS connected to %s. The digital gate is always internal to the Pulseblaster.'%(self.parent_device.name))
            self.gate = DigitalQuantity(self.name+'_gate',self,'gate')
            
    def setamp(self,t,value,units=None):
        self.amplitude.constant(t,value,units)
    def setfreq(self,t,value,units=None):
        self.frequency.constant(t,value,units)
    def setphase(self,t,value,units=None):
        self.phase.constant(t,value,units)
    def enable(self,t):
        if self.gate:
            self.gate.go_high(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.'%(self.name))
                        
    def disable(self,t):
        if self.gate:
            self.gate.go_low(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.'%(self.name))
                
                    
class StaticDDS(Device):
    description = 'Static RF'
    allowed_children = [StaticAnalogQuantity,DigitalOut]
    generation = 2
    def __init__(self,name,parent_device,connection,digital_gate = {},freq_limits = None,freq_conv_class = None,freq_conv_params = {},amp_limits=None,amp_conv_class = None,amp_conv_params = {},phase_limits=None,phase_conv_class = None,phase_conv_params = {}):
        self.clock_type = parent_device.clock_type
        Device.__init__(self,name,parent_device,connection)
        self.frequency = StaticAnalogQuantity(self.name+'_freq',self,'freq',freq_limits,freq_conv_class,freq_conv_params)
        self.amplitude = StaticAnalogQuantity(self.name+'_amp',self,'amp',amp_limits,amp_conv_class,amp_conv_params)
        self.phase = StaticAnalogQuantity(self.name+'_phase',self,'phase',phase_limits,phase_conv_class,phase_conv_params)        
        if isinstance(self.parent_device,NovaTechDDS9M):
            self.frequency.default_value = 0.1
            if 'device' in digital_gate and 'connection' in digital_gate:            
                self.gate = DigitalOut(self.name+'_gate',digital_gate['device'],digital_gate['connection'])
            # Did they only put one key in the dictionary, or use the wrong keywords?
            elif len(digital_gate) > 0:
                raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.'%(self.name))
                
    def setamp(self,value,units=None):
        self.amplitude.constant(value,units)
        
    def setfreq(self,value,units=None):
        self.frequency.constant(value,units)
        
    def setphase(self,value,units=None):
        self.phase.constant(value,units) 
               
    def enable(self,t):
        if self.gate:
            self.gate.go_high(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.'%(self.name))
                        
    def disable(self,t):
        if self.gate:
            self.gate.go_low(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.'%(self.name))
              
            
class RFBlaster(PseudoClock):
    description = 'RF Blaster Rev1.1'
    clock_limit = 500e3
    clock_type = 'fast clock'
    allowed_children = [DDS]
    
    # This will be replaced when we code up arbitrary triggering of RFblasters
    trigger_time = 0.0
    
    def __init__(self,name,ip_address):
        PseudoClock.__init__(self,name,None,None)
        self.BLACS_connection = ip_address
    
    def generate_code(self, hdf5_file):
        from rfblaster import caspr
        import rfblaster.rfjuice
        rfjuice_folder = os.path.dirname(rfblaster.rfjuice.__file__)
        
        import rfblaster.rfjuice.const as c
        from rfblaster.rfjuice.cython.make_diff_table import make_diff_table
        from rfblaster.rfjuice.cython.compile import compileD
        # from rfblaster.rfjuice.compile import compileD
        import tempfile
        from subprocess import Popen, PIPE
        
        # Generate clock and save raw instructions to the h5 file:
        PseudoClock.generate_code(self, hdf5_file)
        dtypes = [('time',float),('amp0',float),('freq0',float),('phase0',float),('amp1',float),('freq1',float),('phase1',float)]
        data = zeros(len(self.times),dtype=dtypes)
        data['time'] = self.times
        for dds in self.child_devices:
            prefix, connection = dds.connection.split()
            data['freq%s'%connection] = dds.frequency.raw_output
            data['amp%s'%connection] = dds.amplitude.raw_output
            data['phase%s'%connection] = dds.phase.raw_output
        group = hdf5_file['devices'].create_group(self.name)
        group.create_dataset('TABLE_DATA',data=data)
        
        # Quantise the data and save it to the h5 file:
        quantised_dtypes = [('time',int32),('amp0',int32),('freq0',int32),('phase0',int32),('amp1',int32),('freq1',int32),('phase1',int32)]
        quantised_data = zeros(len(self.times),dtype=quantised_dtypes)
        quantised_data['time'] = array(c.tT*1e6*data['time']+0.5)
        for dds in range(2):
            # TODO: bounds checking
            # Adding 0.5 to each so that casting to integer rounds:
            quantised_data['freq%d'%dds] = array(c.fF*1e-6*data['freq%d'%dds] + 0.5)
            quantised_data['amp%d'%dds]  = array((2**c.bitsA - 1)*data['amp%d'%dds] + 0.5)
            quantised_data['phase%d'%dds] = array(c.pP*data['amp%d'%dds] + 0.5)
        group.create_dataset('QUANTISED_DATA',data=quantised_data)
        # Generate some assembly code and compile it to machine code:
        assembly_group = group.create_group('ASSEMBLY_CODE')
        binary_group = group.create_group('BINARY_CODE')
        diff_group = group.create_group('DIFF_TABLES')
        for dds in range(2):
            abs_table = zeros((len(self.times), 4),dtype=int32)
            abs_table[:,0] = quantised_data['time']
            abs_table[:,1] = quantised_data['amp%d'%dds]
            abs_table[:,2] = quantised_data['freq%d'%dds]
            abs_table[:,3] = quantised_data['phase%d'%dds]
            
#            abs_table_0 = abs_table.copy()[:-1]
#            abs_table_1 = abs_table.copy()[:-1]
#            abs_table_2 = abs_table.copy()[:-1]
#            abs_tables = [abs_table_0, abs_table_1, abs_table_2]
            
            abs_tables = [abs_table]
            # convert to diff tables:
            diff_tables = [make_diff_table(tab) for tab in abs_tables]
            # Create temporary files, get their paths, and close them:
            with tempfile.NamedTemporaryFile(delete=False) as f:
                temp_assembly_filepath = f.name
            with tempfile.NamedTemporaryFile(delete=False) as f:
                temp_binary_filepath = f.name
                
            try:
                # Compile to assembly:
                with open(temp_assembly_filepath,'w') as assembly_file:
                    for i, dtab in enumerate(diff_tables):
                        compileD(dtab, assembly_file, init=(i == 0),
                                 jump_to_start=(i == 0),
                                 jump_from_end=(i == len(diff_tables) - 1),
                                 local_loop_pre = str(i),
                                 set_defaults = False)
                # Save the assembly to the h5 file:
                with open(temp_assembly_filepath,) as assembly_file:
                    assembly_code = assembly_file.read()
                    assembly_group.create_dataset('DDS%d'%dds, data=assembly_code)
                    for i, diff_table in enumerate(diff_tables):
                        diff_group.create_dataset('DDS%d_difftable%d'%(dds,i), data=diff_table)
                # compile to binary:
                compilation = Popen([caspr,temp_assembly_filepath,temp_binary_filepath],
                                     stdout=PIPE, stderr=PIPE, cwd=rfjuice_folder,startupinfo=startupinfo)
                stdout, stderr = compilation.communicate()
                if compilation.returncode:
                    print stdout
                    raise LabscriptError('RFBlaster compilation exited with code %d\n\n'%compilation.returncode + 
                                         'Stdout was:\n %s\n'%stdout + 'Stderr was:\n%s\n'%stderr)
                # Save the binary to the h5 file:
                with open(temp_binary_filepath,'rb') as binary_file:
                    binary_data = binary_file.read()
                binary_group.create_dataset('DDS%d'%dds, data=binary_data)
            finally:
                # Delete the temporary files:
                os.remove(temp_assembly_filepath)
                os.remove(temp_binary_filepath)
                # print 'assembly:', temp_assembly_filepath
                # print 'binary for dds %d on %s:'%(dds,self.name), temp_binary_filepath
        
                            
class NovaTechDDS9M(IntermediateDevice):
    description = 'NT-DDS9M'
    allowed_children = [DDS, StaticDDS]
    clock_limit = 500e3 # TODO: find out what the actual max clock rate is.
    
    def __init__(self, name, parent_device, clock_type, com_port):
        IntermediateDevice.__init__(self, name, parent_device,clock_type)
        self.BLACS_connection = com_port
    
    def quantise_freq(self,data, device):
        # Ensure that frequencies are within bounds:
        if any(data > 171e6 )  or any(data < 0.1 ):
            raise LabscriptError('%s %s '%(device.description, device.name) +
                              'can only have frequencies between 0.1Hz and 171MHz, ' + 
                              'the limit imposed by %s.'%self.name)
        # It's faster to add 0.5 then typecast than to round to integers first:
        data = array((10*data)+0.5,dtype=uint32)
        scale_factor = 10
        return data, scale_factor
        
    def quantise_phase(self,data,device):
        # ensure that phase wraps around:
        data %= 360
        # It's faster to add 0.5 then typecast than to round to integers first:
        data = array((45.511111111111113*data)+0.5,dtype=uint16)
        scale_factor = 45.511111111111113
        return data, scale_factor
        
    def quantise_amp(self,data,device):
        # ensure that amplitudes are within bounds:
        if any(data > 1 )  or any(data < 0):
            raise LabscriptError('%s %s '%(device.description, device.name) +
                              'can only have amplitudes between 0 and 1 (Volts peak to peak approx), ' + 
                              'the limit imposed by %s.'%self.name)
        # It's faster to add 0.5 then typecast than to round to integers first:
        data = array((1023*data)+0.5,dtype=uint16)
        scale_factor = 1023
        return data, scale_factor
        
    def generate_code(self, hdf5_file):
        DDSs = {}
        for output in self.child_devices:
            # Check that the instructions will fit into RAM:
            if isinstance(output, DDS) and len(output.frequency.raw_output) > 16384 - 2: # -2 to include space for dummy instructions
                raise LabscriptError('%s can only support 16383 instructions. '%self.name +
                                     'Please decrease the sample rates of devices on the same clock, ' + 
                                     'or connect %s to a different pseudoclock.'%self.name)
            try:
                prefix, channel = output.connection.split()
                channel = int(channel)
            except:
                raise LabscriptError('%s %s has invalid connection string: \'%s\'. '%(output.description,output.name,str(output.connection)) + 
                                     'Format must be \'channel n\' with n from 0 to 4.')
            DDSs[channel] = output
        for connection in DDSs:
            if connection in range(4):
                dds = DDSs[connection]   
                dds.frequency.raw_output, dds.frequency.scale_factor = self.quantise_freq(dds.frequency.raw_output, dds)
                dds.phase.raw_output, dds.phase.scale_factor = self.quantise_phase(dds.phase.raw_output, dds)
                dds.amplitude.raw_output, dds.amplitude.scale_factor = self.quantise_amp(dds.amplitude.raw_output, dds)                   
            else:
                raise LabscriptError('%s %s has invalid connection string: \'%s\'. '%(dds.description,dds.name,str(dds.connection)) + 
                                     'Format must be \'channel n\' with n from 0 to 4.')
                                
        dtypes = [('freq%d'%i,uint32) for i in range(2)] + \
                 [('phase%d'%i,uint16) for i in range(2)] + \
                 [('amp%d'%i,uint16) for i in range(2)]
                 
        static_dtypes = [('freq%d'%i,uint32) for i in range(2,4)] + \
                        [('phase%d'%i,uint16) for i in range(2,4)] + \
                        [('amp%d'%i,uint16) for i in range(2,4)]
                        
        if self.clock_type == 'slow clock':
            times = self.parent_device.change_times
        else:
            times = self.parent_device.times
        out_table = zeros(len(times),dtype=dtypes)
        out_table['freq0'].fill(1)
        out_table['freq1'].fill(1)
        
        static_table = zeros(1, dtype=static_dtypes)
        static_table['freq2'].fill(1)
        static_table['freq3'].fill(1)
        
        for connection in range(2):
            if not connection in DDSs:
                continue
            dds = DDSs[connection]
            # The last two instructions are left blank, for BLACS
            # to fill in at program time.
            out_table['freq%d'%connection][:] = dds.frequency.raw_output
            out_table['amp%d'%connection][:] = dds.amplitude.raw_output
            out_table['phase%d'%connection][:] = dds.phase.raw_output
        for connection in range(2,4):
            if not connection in DDSs:
                continue
            dds = DDSs[connection]
            static_table['freq%d'%connection] = dds.frequency.raw_output[0]
            static_table['amp%d'%connection] = dds.amplitude.raw_output[0]
            static_table['phase%d'%connection] = dds.phase.raw_output[0]
            
        grp = hdf5_file.create_group('/devices/'+self.name)
        grp.attrs['frequency_scale_factor'] = 10
        grp.attrs['amplitude_scale_factor'] = 1023
        grp.attrs['phase_scale_factor'] = 45.511111111111113
        grp.create_dataset('TABLE_DATA',compression=config.compression,data=out_table) 
        grp.create_dataset('STATIC_DATA',compression=config.compression,data=static_table) 


class ZaberStageTLSR150D(StaticAnalogQuantity):
    minval=0
    maxval=76346
    description = 'Zaber Stage T-LSR150D'
    
class ZaberStageTLSR300D(StaticAnalogQuantity):
    minval=0
    maxval=151937
    description = 'Zaber Stage T-LSR300D'
    
class ZaberStageTLS28M(StaticAnalogQuantity):
    minval=0
    maxval=282879
    description = 'Zaber Stage T-LS28-M'


class ZaberStageController(Device):
    allowed_children = [ZaberStageTLSR150D,ZaberStageTLSR300D,ZaberStageTLS28M]
    generation = 0
    def __init__(self, name,com_port):
        Device.__init__(self, name, None, None)
        self.clock_type = None
        self.BLACS_connection = com_port
        
    def generate_code(self, hdf5_file):
        data_dict = {}
        for stage in self.child_devices:
            # Call these functions to finalise the stage, they are standard functions of all subclasses of Output:
            ignore = stage.get_change_times()
            stage.make_timeseries([])
            stage.expand_timeseries()
            connection = [int(s) for s in stage.connection.split() if s.isdigit()][0]
            value = stage.raw_output[0]
            if not stage.minval <= value <= stage.maxval:
                # error, out of bounds
                raise LabscriptError('%s %s has value out of bounds. Set value: %s Allowed range: %s to %s.'%(stage.description,stage.name,str(value),str(stage.minval),str(stage.maxval)))
            if not connection > 0:
                # error, invalid connection number
                raise LabscriptError('%s %s has invalid connection number: %s'%(stage.description,stage.name,str(stage.connection)))
            data_dict[str(stage.connection)] = value
        dtypes = [(conn, int) for conn in data_dict]
        data_array = zeros(1, dtype=dtypes)
        for conn in data_dict:
            data_array[0][conn] = data_dict[conn] 
        grp = hdf5_file.create_group('/devices/'+self.name)
        grp.create_dataset('static_values', data=data_array)
        

class LabscriptError(Exception):
    pass
            
def generate_connection_table(hdf5_file):
    connection_table = []
    devicedict = {}
    def sortkey(row):
        device = devicedict[row[0]]
        return str(device.generation) + device.name
    
    # This starts at 4 to accomodate "None"
    max_cal_param_length = 4
    max_BLACS_conn_length = 1
    for device in compiler.inventory:
        devicedict[device.name] = device
        
        # If the device has calibration parameters, then run some checks
        if hasattr(device,"unit_conversion_parameters"):
            try:
                # Are we able to store the calibration parameter dictionary in the h5 file as a string?
                assert(eval(repr(device.unit_conversion_parameters)) == device.unit_conversion_parameters)
            except(AssertionError,SyntaxError):
                raise LabscriptError('The calibration parameters for device "%s" are too complex to store as a string in the connection table'%device.name)
                            
            # Find the logest parameter string
            cal_params = repr(device.unit_conversion_parameters)
            if len(cal_params) > max_cal_param_length:
                max_cal_param_length = len(cal_params)
        else:
            cal_params = str(None)
        
        
        # If the device has a BLACS_connection atribute, then check to see if it is longer than the size of the hdf5 column
        if hasattr(device,"BLACS_connection"):
            # Make sure it is a string!
            BLACS_connection = str(device.BLACS_connection)
            if len(BLACS_connection) > max_BLACS_conn_length:
                max_BLACS_conn_length = len(BLACS_connection)
        else:
            BLACS_connection = ""
            
        connection_table.append((device.name, device.__class__.__name__,
                                 device.parent_device.name if device.parent_device else str(None),
                                 str(device.connection if device.parent_device else str(None)),
                                 device.unit_conversion_class.__name__ if hasattr(device,"unit_conversion_class") and device.unit_conversion_class is not None else str(None),
                                 cal_params,
                                 BLACS_connection))
    
    connection_table.sort(key=sortkey)
    connection_table_dtypes = [('name','a256'), ('class','a256'), ('parent','a256'), ('parent port','a256'),
                               ('unit conversion class','a256'), ('unit conversion params','a'+str(max_cal_param_length)), 
                               ('BLACS_connection','a'+str(max_BLACS_conn_length))]
    connection_table_array = empty(len(connection_table),dtype=connection_table_dtypes)
    for i, row in enumerate(connection_table):
        connection_table_array[i] = row
    hdf5_file.create_dataset('connection table', compression=config.compression, data=connection_table_array, maxshape=(None,))
  
  
def save_labscripts(hdf5_file):
    if compiler.labscript_file is not None:
        script_text = open(compiler.labscript_file).read()
    else:
        script_text = ''
    script = hdf5_file.create_dataset('script',compression=config.compression,data=script_text)
    script.attrs['name'] = os.path.basename(compiler.labscript_file) if compiler.labscript_file is not None else ''
    script.attrs['path'] = os.path.dirname(compiler.labscript_file) if compiler.labscript_file is not None else sys.path[0]
    try:
        import labscriptlib
        prefix = os.path.dirname(labscriptlib.__file__)
        for module in sys.modules.values():
            if hasattr(module,'__file__'):
                path = os.path.abspath(module.__file__)
                if path.startswith(prefix) and (path.endswith('.pyc') or path.endswith('.py')):
                    path = path.replace('.pyc','.py')
                    save_path = 'labscriptlib/' + path.replace(prefix,'').replace('\\','/')
                    if save_path in hdf5_file:
                        # Don't try to save the same module script twice! 
                        # (seems to at least double count __init__.py when you import an entire module as in from labscriptlib.stages import * where stages is a folder with an __init__.py file.
                        # Doesn't seem to want to double count files if you just import the contents of a file within a module
                        continue
                    hdf5_file.create_dataset(save_path, compression=config.compression, data=open(path).read())
                    process = subprocess.Popen(['svn', 'info', path], stdout=subprocess.PIPE,stderr=subprocess.PIPE,startupinfo=startupinfo)
                    info, err = process.communicate()
                    hdf5_file[save_path].attrs['svn info'] = info + '\n' + err
    except ImportError:
        pass
    except WindowsError if os.name == 'nt' else None:
        sys.stderr.write('Warning: Cannot save SVN data for imported scripts. Check that the svn command can be run from the command line')
        
      
def generate_code():
    if compiler.hdf5_filename is None:
        raise LabscriptError('hdf5 file for compilation not set. Please call labscript_init')
    with h5py.File(compiler.hdf5_filename) as hdf5_file:
        hdf5_file.create_group('devices')
        hdf5_file.create_group('calibrations')
        for device in compiler.inventory:
            if not device.parent_device:
                device.generate_code(hdf5_file)
        generate_connection_table(hdf5_file)
        save_labscripts(hdf5_file)
                           
def stop(t):
    # Indicate the end of an experiment and initiate compilation:
    if t == 0:
        raise LabscriptError('Stop time cannot be t=0. Please make your run a finite duration')
    for device in compiler.inventory:
        if not device.parent_device:  
            device.stop_time = t
    generate_code()


def load_globals(hdf5_filename):
    with h5py.File(hdf5_filename,'r') as hdf5_file:
        params = dict(hdf5_file['globals'].attrs)
        for name in params.keys():
            if name in globals().keys() or name in locals().keys() or name in dir(__builtins__):
                raise LabscriptError('Error whilst parsing globals from %s. \'%s\''%(hdf5_filename,name) +
                                     ' is already a name used by Python, labscript, or Pylab.'+
                                     ' Please choose a different variable name to avoid a conflict.')
            if name in keyword.kwlist:
                raise LabscriptError('Error whilst parsing globals from %s. \'%s\''%(hdf5_filename,name) +
                                     ' is a reserved Python keyword.' +
                                     ' Please choose a different variable name.')
            try:
                assert '.' not in name
                exec(name + ' = 0')
                exec('del ' + name )
            except:
                raise LabscriptError('ERROR whilst parsing globals from %s. \'%s\''%(hdf5_filename,name) +
                                     'is not a valid Python variable name.' +
                                     ' Please choose a different variable name.')
                                     
            # Workaround for the fact that numpy.bool_ objects dont 
            # match python's builtin True and False when compared with 'is':
            if type(params[name]) == bool_: # bool_ is numpy.bool_, imported from pylab
                params[name] = bool(params[name])                         
            
            __builtins__[name] = params[name]
            
            
def labscript_init(hdf5_filename, labscript_file=None, new=False):
    # Backup builtins so they can be restored at deinit:
    if new:
        with h5py.File(hdf5_filename ,'w') as hdf5_file:
            hdf5_file.create_group('globals')
    if not os.path.exists(hdf5_filename):
        raise LabscriptError('Provided hdf5 filename %s doesn\'t exist.'%hdf5_filename)
    load_globals(hdf5_filename)
    compiler.hdf5_filename = hdf5_filename
    compiler.labscript_file = os.path.abspath(labscript_file)

def labscript_cleanup():
    """restores builtins and the labscript module to its state before
    labscript_init() was called"""
    for name in _builtins_dict.copy():
        if name not in _existing_builtins_dict:
            del _builtins_dict[name]
        else:
            _builtins_dict[name] = _existing_builtins_dict[name]
    compiler.inventory = []
    compiler.hdf5_filename = None
    compiler.labscript_file = None
    
class compiler:
    # The labscript file being compiled:
    labscript_file = None
    # All defined devices:
    inventory = []
    # The filepath of the h5 file containing globals and which will
    # contain compilation output:
    hdf5_filename = None

   
# Initialisation, runs at import. Can be suppressed by setting
# labscript_auto_init = False in your __main__ module before importing
# labscript. If you do this, you'll need to call labscript_init() yourself:
if not hasattr(__main__, 'labscript_auto_init') or __main__.labscript_auto_init == True:
    if len(sys.argv) > 1:
        labscript_init(sys.argv[1],labscript_file=sys.argv[0])
    elif sys.argv[0]:
        labscript_init(sys.argv[0].replace('.py','.h5'), labscript_file=sys.argv[0], new=True)
