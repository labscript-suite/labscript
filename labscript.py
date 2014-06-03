#####################################################################
#                                                                   #
# /labscript.py                                                     #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

from __future__ import division
import os
import sys
import subprocess
import keyword

import labscript_utils.h5_lock, h5py
from pylab import *

import functions
try:
    from labscript_utils.unitconversions import *
except:
    print 'Failed to import unit conversion classes'

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
# when labscript_cleanup() is called. 
import __builtin__
_builtins_dict = __builtin__.__dict__
_existing_builtins_dict = _builtins_dict.copy()
    
# Startupinfo, for ensuring subprocesses don't launch with a visible command window:
if os.name=='nt':
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= 1 #subprocess.STARTF_USESHOWWINDOW # This variable isn't defined, but apparently it's equal to one.
else:
    startupinfo = None
        
        
class config:
    suppress_mild_warnings = True
    suppress_all_warnings = False
    compression = None # set to 'gzip' for compression 
   
    
class NoWarnings(object):
    """A context manager which sets config.suppress_mild_warnings to True
    whilst in use.  Allows the user to suppress warnings for specific
    lines when they know that the warning does not indicate a problem."""
    def __enter__(self):
        self.existing_warning_setting = config.suppress_all_warnings
        config.suppress_all_warnings = True
    def __exit__(self, *args):
        config.suppress_all_warnings = self.existing_warning_setting
    
no_warnings = NoWarnings() # This is the object that should be used, not the class above

def max_or_zero(*args, **kwargs):
    """returns max(*args) or zero if given an empty sequence (in which case max() would throw an error)"""
    if not args:
        return 0
    if not args[0]:
        return 0
    else:
        return max(*args, **kwargs)
    
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
        if parent_device:
            parent_device.add_device(self)
        
        # Check that the name doesn't already exist in the python namespace
        if name in locals() or name in globals() or name in _builtins_dict:
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
        _builtins_dict[name] = self
        
        # Add self to the compiler's device inventory
        compiler.inventory.append(self)
        
    def add_device(self,device):
        if any([isinstance(device,DeviceClass) for DeviceClass in self.allowed_children]):
            self.child_devices.append(device)
        else:
            raise LabscriptError('Devices of type %s cannot be attached to devices of type %s.'%(device.description,self.description))
    
    @property    
    def pseudoclock(self):
        if isinstance(self, PseudoClock):
            return self
        parent = self.parent_device
        try:
            while not isinstance(parent,PseudoClock):
                parent = parent.parent_device
            return parent
        except Exception as e:
            raise LabscriptError('Couldn\'t find parent pseudoclock of %s, what\'s going on? Original error was %s.'%(self.name, str(e)))
    
    @property
    def t0(self):
        """The earliest time output can be commanded from this device at the start of the experiment.
        This is nonzeo on secondary pseudoclocks due to triggering delays."""
        parent = self.pseudoclock
        if parent.is_master_pseudoclock:
            return 0
        else:
            return round(parent.trigger_times[0] + parent.trigger_delay, 10)
                            
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
    trigger_edge_type = 'rising'
    # How long after a trigger the next instruction is actually output:
    trigger_delay = 0
    # How long a trigger line must remain high/low in order to be detected:
    trigger_minimum_duration = 0 
    # How long after the start of a wait instruction the device is actually capable of resuming:
    wait_delay = 0
    
    def __init__(self,name,trigger_device=None,trigger_connection=None):
        if trigger_device is None:
            parent_device = None
            connection = None
            for device in compiler.inventory:
                if isinstance(device,PseudoClock) and device.is_master_pseudoclock:
                    raise LabscriptError('There is already a master pseudoclock: %s.'%device.name + 
                                         'There cannot be multiple master pseudoclocks - please provide a trigger_device for one of them.')
        else:
            # Instantiate a Trigger object to handle the generation
            # of edges for triggering this clock.  A trigger object is
            # basically just a digital output:
            self.trigger_device = Trigger(name + '_trigger', trigger_device, trigger_connection, self.trigger_edge_type)
            parent_device = self.trigger_device
            connection = 'trigger'
            # Ensure that the parent pseudoclock is, in fact, the master pseudoclock.
            if not trigger_device.pseudoclock.is_master_pseudoclock:
                raise LabscriptError('All secondary pseudoclocks must be triggered by a device being clocked by the master pseudoclock.' +
                                     'Pseudoclocks triggering each other in series is not supported.')
        Device.__init__(self, name, parent_device, connection)
        self.trigger_times = []
        self.wait_times = []
        self.initial_trigger_time = 0
    
    @property    
    def is_master_pseudoclock(self):
        return self.parent_device is None
    
    def set_initial_trigger_time(self, t):
        if compiler.start_called:
            raise LabscriptError('Initial trigger times must be set prior to calling start()')
        if self.is_master_pseudoclock:
            raise LabscriptError('Initial trigger time of master clock is always zero, it cannot be changed.')
        else:
            self.initial_trigger_time = t
            
    def trigger(self, t, duration, wait_delay = 0):
        """Ask the trigger device to produce a digital pulse of a given duration to trigger this pseudoclock"""
        if t == 'initial':
            t = self.initial_trigger_time
        t = round(t,10)
        if self.is_master_pseudoclock:
            if compiler.wait_monitor is not None:
                # Make the wait monitor pulse to signify starting or resumption of the experiment:
                compiler.wait_monitor.trigger(t, duration)
            elif t != self.initial_trigger_time:
                raise LabscriptError("You cannot use waits in unless you have a wait monitor." +
                                     "Please instantiate a WaitMonitor in your connection table.")
            self.trigger_times.append(t)
        else:
            self.trigger_device.trigger(t, duration)
            self.trigger_times.append(round(t + wait_delay,10))
    
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
        if not self.stop_time in change_times:
            dt = self.stop_time - change_times[-1]
            if abs(dt) < 1.0/self.clock_limit:
                raise LabscriptError('The stop time of the experiment is t= %s s, but the last instruction for a device attached to %s is at t= %s s. '%( str(self.stop_time), self.name, str(change_times[-1])) +
                                     'One or more connected devices cannot support update delays shorter than %s sec. Please set the stop_time a bit later.'%str(1.0/self.clock_limit))
            change_times.append(self.stop_time)
        # include trigger times in change_times, so that pseudoclocks always have an instruction immediately following a wait:
        change_times.extend(self.trigger_times)
        # Get rid of duplicates if trigger times were already in the list:
        change_times = list(set(change_times))
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
        all_times = {}
        clocks_in_use = []
        for output in outputs:
            if output.parent_device.clock_type != 'slow clock':            
                if output.parent_device.clock_type not in all_times:
                    all_times[output.parent_device.clock_type] = []
                if output.parent_device.clock_type not in clocks_in_use:
                    clocks_in_use.append(output.parent_device.clock_type)
        
        clock = []
        for i, time in enumerate(change_times):
            if time in self.trigger_times[1:]:
                # A wait instruction:
                clock.append('WAIT')
                
            # list of enabled clocks
            enabled_clocks = []
            
            # what's the fastest clock rate?
            maxrate = 0
            for output in outputs:
                # Check if output is sweeping and has highest clock rate
                # so far. If so, store its clock rate to max_rate:
                if hasattr(output,'timeseries') and isinstance(output.timeseries[i],dict):
                    # find out what clock this is on so we can enable it
                    enabled_clocks.append(output.parent_device.clock_type)
                    if output.timeseries[i]['clock rate'] > maxrate:
                        # It does have the highest clock rate? Then store that rate to max_rate:
                        maxrate = output.timeseries[i]['clock rate']
            if maxrate:
                # round to the nearest clock rate that the pseudoclock can actually support:
                period = 1/maxrate
                quantised_period = period/self.clock_resolution
                quantised_period = round(quantised_period)
                period = quantised_period*self.clock_resolution
                maxrate = 1/period
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
                
                for a_clock in clocks_in_use:
                    if a_clock in enabled_clocks:  
                        all_times[a_clock].append(ticks)
                    else:
                        all_times[a_clock].append(time)
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
                        clock.append({'start': time, 'reps': 1, 'step': 1/float(maxrate),'slow_clock_tick':True, 'fast_clock':enabled_clocks})
                        clock.append({'start': time + 1/float(maxrate), 'reps': n_ticks-2, 'step': 1/float(maxrate),'slow_clock_tick':False, 'fast_clock':enabled_clocks})
                    else:
                        clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate),'slow_clock_tick':True, 'fast_clock':enabled_clocks})
                # The last clock tick has a different duration depending
                # on the next step. The slow clock must tick here if it
                # hasn't ticked already, that is if n_ticks = 1.
                clock.append({'start': ticks[-1], 'reps': 1, 'step': change_times[i+1] - ticks[-1],'slow_clock_tick': True if n_ticks == 1 else False, 'fast_clock':enabled_clocks})
            else:
                for a_clock in clocks_in_use:
                    all_times[a_clock].append(time)
                try: 
                    # If there was no ramping, here is a single clock tick:
                    clock.append({'start': time, 'reps': 1, 'step': change_times[i+1] - time,'slow_clock_tick':True, 'fast_clock':'all'})
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
                        clock.append({'start': time, 'reps': 1, 'step': 10.0/self.clock_limit,'slow_clock_tick':True, 'fast_clock':'all'})
        return all_times, clock
    
    def do_checks(self, outputs):
        """Basic error checking te ensure the user's instructions make sense"""
        for output in outputs:
            output.do_checks(self.trigger_times)
            
    def offset_instructions_from_trigger(self, outputs):
        for output in outputs:
            output.offset_instructions_from_trigger(self.trigger_times)
        
        if not self.is_master_pseudoclock:
            # Adjust the stop time relative to the last trigger time
            self.stop_time = self.stop_time - self.trigger_delay * len(self.trigger_times)
            # Modify the trigger times themselves so that we insert wait instructions at the right times:
            initial_trigger_time = self.trigger_times[0]
            self.trigger_times = [t - initial_trigger_time - i*self.trigger_delay for i, t in enumerate(self.trigger_times)]
                            
    def generate_clock(self):
        outputs = self.get_all_outputs()
        self.do_checks(outputs)
        self.offset_instructions_from_trigger(outputs)
        change_times = self.collect_change_times(outputs)
        for output in outputs:
            output.make_timeseries(change_times)
        all_times, clock = self.expand_change_times(change_times, outputs)
        for output in outputs:
            if output.parent_device.clock_type == 'slow clock':
                output.expand_timeseries(change_times)
            else:
                output.expand_timeseries(all_times[output.parent_device.clock_type])
        self.clock = clock
        self.change_times = fastflatten(change_times, float)
        
        # Dummy empty arrays in case there are no outputs on one or both clocks:
        self.times = {}
        for clock_type, time_array in all_times.items():
            self.times[clock_type] = fastflatten(time_array,float)
        
    def generate_code(self, hdf5_file):
        self.generate_clock()
        Device.generate_code(self, hdf5_file)
        
        
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
    
    @property
    def clock_limit(self):
        parent = self.pseudoclock
        return parent.clock_limit
    
    @property
    def trigger_delay(self):
        """The earliest time output can be commanded from this device after a trigger.
        This is nonzeo on secondary pseudoclocks due to triggering delays."""
        parent = self.pseudoclock
        if parent.is_master_pseudoclock:
            return 0
        else:
            return parent.trigger_delay
    
    @property
    def wait_delay(self):
        """The earliest time output can be commanded from this device after a wait.
        This is nonzeo on secondary pseudoclocks due to triggering delays and the fact
        that the master clock doesn't provide a resume trigger to secondary clocks until
        a minimum time has elapsed: compiler.wait_delay. This is so that if a wait is 
        extremely short, the child clock is actually ready for the trigger.
        """
        delay = compiler.wait_delay if self.pseudoclock.is_master_pseudoclock else 0
        return self.trigger_delay + delay
            
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
        if not compiler.start_called:
            raise LabscriptError('Cannot add instructions prior to calling start()')
        # round to the nearest 0.1 nanoseconds, to prevent floating point
        # rounding errors from breaking our equality checks later on.
        time = round(time,10)
        # Also round end time of ramps to the nearest 0.1 ns:
        if isinstance(instruction,dict):
            instruction['end time'] = round(instruction['end time'],10)
            instruction['initial time'] = round(instruction['initial time'],10)
        # Check that time is not negative or too soon after t=0:
        if time < self.t0:
            err = ' '.join([self.description, self.name, 'has an instruction at t=%ss,'%str(time),
                 'Due to the delay in triggering its pseudoclock, the earliest output possible is at t=%s.'%str(self.t0)])
            raise LabscriptError(err)
        # Check that this doesn't collide with previous instructions:
        if time in self.instructions.keys():
            if not config.suppress_all_warnings:
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
    
    def do_checks(self, trigger_times):
        """Basic error checking te ensure the user's instructions make sense"""
        # Check if there are no instructions. Generate a warning and insert an
        # instruction telling the output to remain at its default value.
        if not self.instructions:
            if not config.suppress_mild_warnings and not config.suppress_all_warnings:
                sys.stderr.write(' '.join(['WARNING:', self.name, 'has no instructions. It will be set to %s for all time.\n'%self.instruction_to_string(self.default_value)]))
            self.add_instruction(self.t0, self.default_value)  
        # Check if there are no instructions at the initial time. Generate a warning and insert an
        # instruction telling the output to start at its default value.
        if self.t0 not in self.instructions.keys():
            if not config.suppress_mild_warnings and not config.suppress_all_warnings:
               sys.stderr.write(' '.join(['WARNING:', self.name, 'has no initial instruction. It will initially be set to %s.\n'%self.instruction_to_string(self.default_value)]))
            self.add_instruction(self.t0, self.default_value) 
        # Check that ramps have instructions following them.
        # If they don't, insert an instruction telling them to hold their final value.
        for instruction in self.instructions.values():
            if isinstance(instruction, dict) and instruction['end time'] not in self.instructions.keys():
                self.add_instruction(instruction['end time'], instruction['function'](instruction['end time']-instruction['initial time']), instruction['units'])
        # Checks for trigger times:
        for trigger_time in trigger_times:
            for t, instruction in self.instructions.items():
                # Check no ramps are happening at the trigger time:
                if isinstance(instruction, dict) and instruction['initial time'] < trigger_time and instruction['end time'] > trigger_time:
                    err = (' %s %s has a ramp %s from t = %s to %s. ' % (self.description, 
                            self.name, instruction['description'], str(instruction['initial time']), str(instruction['end time'])) +
                           'This overlaps with a trigger at t=%s, and so cannot be performed.' % str(trigger_time))
                    raise LabscriptError(err)
                # Check that nothing is happening during the delay time after the trigger:
                if round(trigger_time,10) < round(t,10) < round(trigger_time + self.trigger_delay, 10):
                    err = (' %s %s has an instruction at t = %s. ' % (self.description, self.name, str(t)) + 
                           'This is too soon after a trigger at t=%s, '%str(trigger_time) + 
                           'the earliest output possible after this trigger is at t=%s'%str(trigger_time + self.trigger_delay))
                    raise LabscriptError(err)
                # Check that there are no instructions too soon before the trigger:
                if 0 < trigger_time - t < max(self.clock_limit, compiler.wait_delay):
                    err = (' %s %s has an instruction at t = %s. ' % (self.description, self.name, str(t)) + 
                           'This is too soon before a trigger at t=%s, '%str(trigger_time) + 
                           'the latest output possible before this trigger is at t=%s'%str(trigger_time - max(self.clock_limit, compiler.wait_delay)))
                           
    def offset_instructions_from_trigger(self, trigger_times):
        """Subtracts self.trigger_delay from all instructions at or after each trigger_time"""
        offset_instructions = {}
        for t, instruction in self.instructions.items():
            # How much of a delay is there for this instruction? That depends how many triggers there are prior to it:
            n_triggers_prior = len([time for time in trigger_times if time < t])
            # The cumulative offset at this point in time:
            offset = self.trigger_delay * n_triggers_prior + trigger_times[0]
            if isinstance(instruction,dict):
                offset_instruction = instruction.copy()
                offset_instruction['end time'] = instruction['end time'] - offset
                offset_instruction['initial time'] = instruction['initial time'] - offset
            else:
                offset_instruction = instruction
                
            offset_instructions[t - offset] = offset_instruction
        self.instructions = offset_instructions
                   
    def get_change_times(self):
        """If this function is being called, it means that the parent
        Pseudoclock has requested a list of times that this output changes
        state."""        
        times = self.instructions.keys()
        times.sort()
        self.times = times
        return times
    
    def make_timeseries(self, change_times):
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
                        midpoints = zeros(1)
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
        
    def sine4_ramp(self,t,duration,initial,final,samplerate,units=None):
        self.add_instruction(t, {'function': functions.sine4_ramp(duration,initial,final), 'description':'sinusoidal ramp',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})   
                
        return duration
        
    def sine4_reverse_ramp(self,t,duration,initial,final,samplerate,units=None):
        self.add_instruction(t, {'function': functions.sine4_reverse_ramp(duration,initial,final), 'description':'sinusoidal ramp',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})   
                
        return duration
    
    def exp_ramp(self, t, duration, initial, final, samplerate, zero=0, trunc=False, trunc_type='linear', units=None):
        if trunc is not False:
            if trunc_type == 'linear':
                trunc_duration = duration*log((initial-zero)/(trunc-zero))/log((initial-zero)/(final-zero))
            elif trunc_type == 'exponential':
                trunc_duration = trunc * duration
                # final = functions.exp_ramp(0, duration, initial, final, zero)(trunc_duration)
            else:
                raise LabscriptError('Truncation type for exp_ramp not supported. Must be either linear or exponential.')
        else:
            trunc_duration = duration
        self.add_instruction(t, {'function': functions.exp_ramp(duration,initial,final,zero), 'description':'exponential ramp',
                             'initial time':t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
        
        return trunc_duration
     
    def exp_ramp_t(self, t, duration, initial, final, time_constant, samplerate, trunc=False, trunc_type='linear', units=None):
        # Exponential ramp set by the time constant. No truncation yet!
        zero = (final-initial*exp(-duration/time_constant)) / (1-exp(-duration/time_constant))
        if trunc is not False:
            if trunc_type == 'linear':
                trunc_duration = time_constant * log((initial-zero)/(trunc-zero))
            elif trunc_type == 'exponential':
                trunc_duration = trunc * duration
            else:
                raise LabscriptError('Truncation type for exp_ramp_t not supported. Must be either linear or exponential.')
        else:
            trunc_duration = duration
        self.add_instruction(t, {'function': functions.exp_ramp_t(duration, initial, final, time_constant), 'description':'exponential ramp with time consntant',
                             'initial time':t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
                
        return trunc_duration
    

    def piecewise_accel_ramp(self,t,duration,initial,final,samplerate, units=None):
        self.add_instruction(t, {'function': functions.piecewise_accel(duration,initial,final), 'description':'piecewise linear accelleration ramp',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': units})   
                
        return duration
    
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
    
    # Redefine __init__ so that you cannot define a limit or calibration for DO
    def __init__(self,name,parent_device,connection):
        Output.__init__(self,name,parent_device,connection)
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 
    
    '''
    This function only works if the DigitalQuantity is on a fast clock
    
    The pulse_sequence parameter should be specified as a list of tuples. 
    Each tuple should be of the form (time,state)
    
    The period parmeter should, in general, be longer than the entire pulse sequence, 
    and defines how long the final tuple should be held for before repeating the pulse sequence.
    
    The pulse sequence specified will be repeated from time t until t+duration.
    
    The samplerate parameter specifies how often to update the output
    
    Note 1: The samplerate should be significantly faster than the smallest time difference between 
    two states in the pulse sequence, or else points in your pulse sequence may never be evaluated.
    
    Note 2: The time points your pulse sequence is evaluated at may be different than you expect,
    if another output changes state between t and t+duration. As such, you should set the samplerate
    high enough that even if this rounding of tie points occurs (to fit in the update required to change the other output)
    your pulse sequence will not be significantly altered)
    '''
    def repeat_pulse_sequence(self,t,duration,pulse_sequence,period,samplerate):
        self.add_instruction(t, {'function': functions.pulse_sequence(pulse_sequence,period), 'description':'pulse sequence',
                                 'initial time':t, 'end time': t + duration, 'clock rate': samplerate, 'units': None})
        
        return duration

class DigitalOut(DigitalQuantity):
    description = 'digital output'

class StaticDigitalQuantity(DigitalQuantity):
    description = 'static digital quantity'
    value_set = False
    
    def go_high(self):
        if not self.value_set:
            self.add_instruction(0,1)
            self.value_set = True
            self.static_value = 1
    def go_low(self):
        if not self.value_set:
            self.add_instruction(0,0) 
            self.value_set = True
            self.static_value = 0
        else:
            raise LabscriptError('%s %s has already been set to %s. It cannot also be set to %s.'%(self.description, self.name, self.instruction_to_string[self.static_value], self.instruction_to_string[value]))
    
    def get_change_times(self):
        if not hasattr(self,'static_value'):
            self.static_value = self.default_value
            
        return [] # Return an empty list as the calling function at the pseudoclock level expects a list
    
    def make_timeseries(self,change_times):
        pass
    
    def expand_timeseries(self,*args,**kwargs):
        self.raw_output = array([self.static_value])

class StaticDigitalOut(StaticDigitalQuantity):
    description = 'static digital output'
        
class AnalogIn(Device):
    description = 'Analog Input'
    generation = 3
    def __init__(self,name,parent_device,connection,scale_factor=1.0,units='Volts'):
         self.acquisitions = []
         self.scale_factor = scale_factor
         self.units=units
         Device.__init__(self,name,parent_device,connection)
   
    def acquire(self,label,start_time,end_time,wait_label='',scale_factor=None,units=None):
        if scale_factor is None:
            scale_factor = self.scale_factor
        if units is None:
            units = self.units
        self.acquisitions.append({'start_time': start_time, 'end_time': end_time,
                                 'label': label, 'wait_label':wait_label, 'scale_factor':scale_factor,'units':units})
        return end_time - start_time
  
  
class IntermediateDevice(Device):
    generation = 1
    
    def __init__(self, name, parent_device,clock_type):
        self.name = name
        clock_valid = False
        if hasattr(parent_device, 'extra_clocks') :
            if clock_type == 'fast clock':
                clock_type = parent_device.extra_clocks[0]
            if clock_type in parent_device.extra_clocks:
                clock_valid = True
        if not clock_type in ['fast clock', 'slow clock'] and not clock_valid:
            raise LabscriptError("Clock type for %s %s can only be 'slow clock' or 'fast clock'."%(self.name,self.description))
        Device.__init__(self, name, parent_device, clock_type)
        self.clock_type = clock_type
        self.parent_device.clock_limit = min([self.parent_device.clock_limit,self.clock_limit])
    
        
class Shutter(DigitalOut):
    description = 'shutter'
    
    def __init__(self,name,parent_device,connection,delay=(0,0),open_state=1):
        DigitalOut.__init__(self, name, parent_device, connection)
        self.open_delay, self.close_delay = delay
        self.open_state = open_state
        if self.open_state == 1:
            self.allowed_states = {0: 'closed', 1: 'open'}
        elif self.open_state == 0:
            self.allowed_states = {1: 'closed', 0: 'open'}
        else:
            raise LabscriptError("Shutter %s wasn't instantiated with open_state = 0 or 1." % self.name)

    # If a shutter is asked to do something at t=0, it cannot start moving
    # earlier than that.  So initial shutter states will have imprecise
    # timing. Not throwing a warning here because if I did, every run
    # would throw a warning for every shutter. The documentation will
    # have to make a point of this.
    def open(self, t):
        if self.open_state == 1:
            self.go_high(t-self.open_delay if t >= self.open_delay else 0)
        elif self.open_state == 0:
            self.go_low(t-self.open_delay if t >= self.open_delay else 0)

    def close(self, t):
        if self.open_state == 1:
            self.go_low(t-self.close_delay if t >= self.close_delay else 0)  
        elif self.open_state == 0:
            self.go_high(t-self.close_delay if t >= self.close_delay else 0)
    
    def generate_code(self, hdf5_file):
        classname = self.__class__.__name__
        calibration_table_dtypes = [('name','a256'), ('open_delay',float), ('close_delay',float)]
        if classname not in hdf5_file['calibrations']:
            hdf5_file['calibrations'].create_dataset(classname, (0,), dtype=calibration_table_dtypes, maxshape=(None,))
        metadata = (self.name,self.open_delay,self.close_delay)
        dataset = hdf5_file['calibrations'][classname]
        dataset.resize((len(dataset)+1,))
        dataset[len(dataset)-1] = metadata
        
        
class Trigger(DigitalOut):
    description = 'trigger device'
    allowed_states = {1:'high', 0:'low'}
    allowed_children = [PseudoClock]
    def __init__(self, name, parent_device, connection, trigger_edge_type='rising'):
        DigitalOut.__init__(self,name,parent_device,connection)
        self.trigger_edge_type = trigger_edge_type
        if self.trigger_edge_type == 'rising':
            self.enable = self.go_high
            self.disable = self.go_low
            self.allowed_states = {1:'enabled', 0:'disabled'}
        elif self.trigger_edge_type == 'falling':
            self.enable = self.go_low
            self.disable = self.go_high
            self.allowed_states = {1:'disabled', 0:'enabled'}
        else:
            raise ValueError('trigger_edge_type must be \'rising\' or \'falling\', not \'%s\'.'%trigger_edge_type)

    def trigger(self, t, duration):
        if t != self.t0 and self.t0 not in self.instructions:
            self.disable(self.t0)
        self.enable(t)
        self.disable(t + duration)


class WaitMonitor(Trigger):
     def __init__(self, name, parent_device, connection, acquisition_device, acquisition_connection, timeout_device, timeout_connection):
        if compiler.wait_monitor is not None:
            raise LabscriptError("Cannot instantiate a second WaitMonitor: there can be only be one in the experiment")
        compiler.wait_monitor = self
        Trigger.__init__(self, name, parent_device, connection, trigger_edge_type='rising')
        if not parent_device.pseudoclock.is_master_pseudoclock:
            raise LabscriptError('The output device for monitoring wait durations must be clocked by the master pseudoclock')
        self.acquisition_device = acquisition_device
        self.acquisition_connection = acquisition_connection 
        self.timeout_device = timeout_device
        self.timeout_connection = timeout_connection 
        
        
class DDS(Device):
    description = 'DDS'
    allowed_children = [AnalogQuantity,DigitalOut,DigitalQuantity] # Adds its own children when initialised
    generation = 2
    def __init__(self, name, parent_device, connection, digital_gate={}, freq_limits=None, freq_conv_class=None, freq_conv_params={},
                 amp_limits=None, amp_conv_class=None, amp_conv_params={}, phase_limits=None, phase_conv_class=None, phase_conv_params = {}):
        self.clock_type = parent_device.clock_type
        Device.__init__(self, name, parent_device, connection)
        
        if isinstance(self.parent_device, NovaTechDDS9M):
            # create some default unit converstion classes if none are specified
            if freq_conv_class is None:
                if 'NovaTechDDS9mFreqConversion' in globals():
                    freq_conv_class = NovaTechDDS9mFreqConversion
            if amp_conv_class is None:
                if 'NovaTechDDS9mAmpConversion' in globals():
                    amp_conv_class = NovaTechDDS9mAmpConversion
        
        self.frequency = AnalogQuantity(self.name + '_freq', self, 'freq', freq_limits, freq_conv_class, freq_conv_params)
        self.amplitude = AnalogQuantity(self.name + '_amp', self, 'amp', amp_limits, amp_conv_class, amp_conv_params)
        self.phase = AnalogQuantity(self.name + '_phase', self, 'phase', phase_limits, phase_conv_class, phase_conv_params)
        self.gate = None
        if isinstance(self.parent_device, NovaTechDDS9M):
            self.frequency.default_value = 0.1
            if 'device' in digital_gate and 'connection' in digital_gate:            
                self.gate = DigitalOut(self.name + '_gate', digital_gate['device'], digital_gate['connection'])
            # Did they only put one key in the dictionary, or use the wrong keywords?
            elif len(digital_gate) > 0:
                raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.' % (self.name))
        elif isinstance(self.parent_device, PulseBlaster):
            if 'device' in digital_gate and 'connection' in digital_gate: 
                raise LabscriptError('You cannot specify a digital gate for a DDS connected to %s. The digital gate is always internal to the Pulseblaster.' % (self.parent_device.name))
            self.gate = DigitalQuantity(self.name + '_gate', self, 'gate')
            
            
    def setamp(self, t, value, units=None):
        self.amplitude.constant(t, value, units)
        
    def setfreq(self, t, value, units=None):
        self.frequency.constant(t, value, units)
        
    def setphase(self, t, value, units=None):
        self.phase.constant(t, value, units)
        
    def enable(self, t):
        if self.gate:
            self.gate.go_high(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.' % (self.name))
            
    def disable(self, t):
        if self.gate:
            self.gate.go_low(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.' % (self.name))
            
    def pulse(self, duration, amplitude, frequency, phase=None, print_summary=True):
        if print_summary:
            print_time(t, '%s pulse at %.4f MHz for %.3f ms' % (self.name, frequency/MHz, duration/ms))
        self.setamp(t, amplitude)
        if not frequency == None:
            self.setfreq(t, frequency)
        if not phase == None:
            self.setphase(t, phase)
        if amplitude != 0:
            self.enable(t)
        self.disable(t)
        self.setamp(t, 0)
        return duration

        
class StaticDDS(Device):
    description = 'Static RF'
    allowed_children = [StaticAnalogQuantity,DigitalOut,StaticDigitalOut]
    generation = 2
    def __init__(self,name,parent_device,connection,digital_gate = {},freq_limits = None,freq_conv_class = None,freq_conv_params = {},amp_limits=None,amp_conv_class = None,amp_conv_params = {},phase_limits=None,phase_conv_class = None,phase_conv_params = {}):
        self.clock_type = parent_device.clock_type
        Device.__init__(self,name,parent_device,connection)
        
        self.frequency = StaticAnalogQuantity(self.name+'_freq',self,'freq',freq_limits,freq_conv_class,freq_conv_params)
        
        if not isinstance(self.parent_device,PhaseMatrixQuickSyn):
            self.amplitude = StaticAnalogQuantity(self.name+'_amp',self,'amp',amp_limits,amp_conv_class,amp_conv_params)
            self.phase = StaticAnalogQuantity(self.name+'_phase',self,'phase',phase_limits,phase_conv_class,phase_conv_params)        
        
        
        
        if isinstance(self.parent_device,NovaTechDDS9M):
            self.frequency.default_value = 0.1
            if 'device' in digital_gate and 'connection' in digital_gate:            
                self.gate = DigitalOut(self.name+'_gate',digital_gate['device'],digital_gate['connection'])
            # Did they only put one key in the dictionary, or use the wrong keywords?
            elif len(digital_gate) > 0:
                raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.'%(self.name))
        
        elif isinstance(self.parent_device,PhaseMatrixQuickSyn):
            self.frequency.default_value = 0.5e9
            self.gate = StaticDigitalOut(self.name+'_gate',self,'gate')
            
            
    def setamp(self,value,units=None):
        if isinstance(self.parent_device,PhaseMatrixQuickSyn):
            raise LabscriptError('You cannot set the amplitude of the QuickSyn')
        else:
            self.amplitude.constant(value,units)
        
    def setfreq(self,value,units=None):
        self.frequency.constant(value,units)
        
    def setphase(self,value,units=None):
        if isinstance(self.parent_device,PhaseMatrixQuickSyn):
            raise LabscriptError('You cannot set the phase of the QuickSyn')
        else:
            self.phase.constant(value,units) 
            
    def enable(self,t=None):        
        if isinstance(self.parent_device,PhaseMatrixQuickSyn):
            self.gate.go_high()
        
        else:
            if self.gate:
                self.gate.go_high(t)
            else:
                raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.'%(self.name))
                            
    def disable(self,t=None):
        if isinstance(self.parent_device,PhaseMatrixQuickSyn):
            self.gate.go_low()
        
        else:
            if self.gate:
                self.gate.go_low(t)
            else:
                raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.'%(self.name))
              
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
    dataset = hdf5_file.create_dataset('connection table', compression=config.compression, data=connection_table_array, maxshape=(None,))
    dataset.attrs['master_pseudoclock'] = compiler.master_pseudoclock.name
  
def save_labscripts(hdf5_file):
    if compiler.labscript_file is not None:
        script_text = open(compiler.labscript_file).read()
    else:
        script_text = ''
    script = hdf5_file.create_dataset('script',data=script_text)
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
                    hdf5_file.create_dataset(save_path, data=open(path).read())
                    process = subprocess.Popen(['svn', 'info', path], stdout=subprocess.PIPE,stderr=subprocess.PIPE,startupinfo=startupinfo)
                    info, err = process.communicate()
                    hdf5_file[save_path].attrs['svn info'] = info + '\n' + err
    except ImportError:
        pass
    except WindowsError if os.name == 'nt' else None:
        sys.stderr.write('Warning: Cannot save SVN data for imported scripts. Check that the svn command can be run from the command line\n')
        
def generate_wait_table(hdf5_file):
    dtypes = [('label','a256'), ('time', float), ('timeout', float)]
    data_array = zeros(len(compiler.wait_table), dtype=dtypes)
    for i, t in enumerate(sorted(compiler.wait_table)):
        label, timeout = compiler.wait_table[t]
        data_array[i] = label, t, timeout
    dataset = hdf5_file.create_dataset('waits', data = data_array)
    if compiler.wait_monitor is not None:
        acquisition_device = compiler.wait_monitor.acquisition_device.name 
        acquisition_connection = compiler.wait_monitor.acquisition_connection
        timeout_device = compiler.wait_monitor.timeout_device.name 
        timeout_connection = compiler.wait_monitor.timeout_connection
    else:
        acquisition_device, acquisition_connection, timeout_device, timeout_connection = '','','',''
    dataset.attrs['wait_monitor_acquisition_device'] = acquisition_device
    dataset.attrs['wait_monitor_acquisition_connection'] = acquisition_connection
    dataset.attrs['wait_monitor_timeout_device'] = timeout_device
    dataset.attrs['wait_monitor_timeout_connection'] = timeout_connection
    
def generate_code():
    if compiler.hdf5_filename is None:
        raise LabscriptError('hdf5 file for compilation not set. Please call labscript_init')
    elif not os.path.exists(compiler.hdf5_filename):
        with h5py.File(compiler.hdf5_filename ,'w') as hdf5_file:
            hdf5_file.create_group('globals')
    with h5py.File(compiler.hdf5_filename) as hdf5_file:
        try:
            hdf5_file.create_group('devices')
            hdf5_file.create_group('calibrations')
        except ValueError:
            # Group(s) already exist - this is not a fresh h5 file, we cannot compile with it:
            raise ValueError('The HDF5 file %s already contains compilation data '%compiler.hdf5_filename +
                             '(possibly partial if due to failed compilation). ' +
                             'Please use a fresh shot file. ' +
                             'If this one was autogenerated by a previous compilation, ' +
                             'and you wish to have a new one autogenerated, '+
                             'simply delete it and run again, or add the \'-f\' command line ' +
                             'argument to automatically overwrite.')
        for device in compiler.inventory:
            if not device.parent_device:
                device.generate_code(hdf5_file)
        generate_connection_table(hdf5_file)
        generate_wait_table(hdf5_file)
        save_labscripts(hdf5_file)

def trigger_all_pseudoclocks(t='initial'):
    # Must wait this long before providing a trigger, in case child clocks aren't ready yet:
    wait_delay = compiler.wait_delay
    if t == 'initial':
        # But not at the start of the experiment:
        wait_delay = 0
    # Trigger them all:
    for pseudoclock in compiler.all_pseudoclocks:
        pseudoclock.trigger(t, compiler.trigger_duration)
    # How long until all devices can take instructions again? The user
    # can command output from devices on the master clock immediately,
    # but unless things are time critical, they can wait this long and
    # know for sure all devices can receive instructions:
    max_delay_time = max_or_zero([pseudoclock.trigger_delay for pseudoclock in compiler.all_pseudoclocks if not pseudoclock.is_master_pseudoclock])
    # On the other hand, perhaps the trigger duration and clock limit of the master clock is
    # limiting when we can next give devices instructions:
    max_delay = max(compiler.trigger_duration + 1.0/compiler.master_pseudoclock.clock_limit, max_delay_time)
    return max_delay + wait_delay
    
def wait(label, t, timeout=5):
    if not str(label):
        raise LabscriptError('Wait must have a name')
    max_delay = trigger_all_pseudoclocks(t)
    if t in compiler.wait_table:
        raise LabscriptError('There is already a wait at t=%s'%str(t))
    if any([label==existing_label for existing_label, timeout in compiler.wait_table.values()]):
        raise LabscriptError('There is already a wait named %s'%str(label))
    compiler.wait_table[t] = str(label), float(timeout)
    return max_delay

def start():
    compiler.start_called = True
    # Get and save some timing info about the pseudoclocks:
    all_pseudoclocks = [device for device in compiler.inventory if isinstance(device, PseudoClock)]
    compiler.all_pseudoclocks = all_pseudoclocks
    master_pseudoclock, = [pseudoclock for pseudoclock in all_pseudoclocks if pseudoclock.is_master_pseudoclock]
    compiler.master_pseudoclock = master_pseudoclock
    # Which pseudoclock requires the longest pulse in order to trigger it?
    compiler.trigger_duration = max_or_zero([pseudoclock.trigger_minimum_duration for pseudoclock in all_pseudoclocks if not pseudoclock.is_master_pseudoclock])
    # Provide this, or the minimum possible pulse, whichever is longer:
    compiler.trigger_duration = max(2.0/master_pseudoclock.clock_limit, compiler.trigger_duration)
    # Must wait this long before providing a trigger, in case child clocks aren't ready yet:
    compiler.wait_delay = max_or_zero([pseudoclock.wait_delay for pseudoclock in all_pseudoclocks if not pseudoclock.is_master_pseudoclock])
    
    # Have the master clock trigger pseudoclocks at t = 0:
    max_delay = trigger_all_pseudoclocks()
    return max_delay
    
def stop(t):
    # Indicate the end of an experiment and initiate compilation:
    if t == 0:
        raise LabscriptError('Stop time cannot be t=0. Please make your run a finite duration')
    for device in compiler.inventory:
        if isinstance(device, PseudoClock):
            device.stop_time = t
    generate_code()

def load_globals(hdf5_filename):
    with h5py.File(hdf5_filename,'r') as hdf5_file:
        params = dict(hdf5_file['globals'].attrs)
        for name in params.keys():
            if name in globals() or name in locals() or name in _builtins_dict:
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
            
            _builtins_dict[name] = params[name]
            
            
def labscript_init(hdf5_filename, labscript_file=None, new=False, overwrite=False):
    if new:
        # defer file creation until generate_code(), so that filesystem
        # is not littered with h5 files when the user merely imports
        # labscript. If the file already exists, and overwrite is true, delete it so we get one fresh.
        if os.path.exists(hdf5_filename) and overwrite:
            os.unlink(hdf5_filename)
    elif not os.path.exists(hdf5_filename):
        raise LabscriptError('Provided hdf5 filename %s doesn\'t exist.'%hdf5_filename)
    else:
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
    compiler.start_called = False
    compiler.wait_table = {}
    compiler.wait_monitor = None
    compiler.master_pseudoclock = None
    compiler.all_pseudoclocks = None
    compiler.trigger_duration = 0
    compiler.wait_delay = 0
    
class compiler:
    # The labscript file being compiled:
    labscript_file = None
    # All defined devices:
    inventory = []
    # The filepath of the h5 file containing globals and which will
    # contain compilation output:
    hdf5_filename = None
    start_called = False
    wait_table = {}
    wait_monitor = None
    master_pseudoclock = None
    all_pseudoclocks = None
    trigger_duration = 0
    wait_delay = 0

