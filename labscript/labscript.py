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

import builtins
import os
import sys
import subprocess
import keyword
from inspect import getcallargs
from functools import wraps

# Notes for v3
#
# Anything commented with TO_DELETE:runmanager-batchompiler-agnostic 
# can be removed with a major version bump of labscript (aka v3+)
# We leave it in now to maintain backwards compatibility between new labscript
# and old runmanager.
# The code to be removed relates to the move of the globals loading code from
# labscript to runmanager batch compiler.


import labscript_utils.h5_lock, h5py
import labscript_utils.properties

# This imports the default Qt library that other labscript suite code will
# import as well, since it all uses qtutils. By having a Qt library already
# imported, we ensure matplotlib (imported by pylab) will notice this and use
# the same Qt library and API version, and therefore not conflict with any
# other code is using:
import qtutils

from labscript_utils import dedent
from labscript_utils.properties import set_attributes

from pylab import *

import labscript.functions as functions
try:
    from labscript_utils.unitconversions import *
except ImportError:
    sys.stderr.write('Warning: Failed to import unit conversion classes\n')

ns = 1e-9
us = 1e-6
ms = 1e-3
s = 1
Hz = 1
kHz = 1e3
MHz = 1e6
GHz = 1e9

# Create a reference to the builtins dict 
# update this if accessing builtins ever changes
_builtins_dict = builtins.__dict__
    
# Startupinfo, for ensuring subprocesses don't launch with a visible command window:
if os.name=='nt':
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= 1 #subprocess.STARTF_USESHOWWINDOW # This variable isn't defined, but apparently it's equal to one.
else:
    startupinfo = None
        
        
class config(object):
    suppress_mild_warnings = True
    suppress_all_warnings = False
    compression = 'gzip'  # set to 'gzip' for compression 
   
    
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
    if np.array_equal(arrays[0], 0):
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

def set_passed_properties(property_names = {}):
    """
    Decorator for device __init__ methods that saves the listed arguments/keyword
    arguments as properties. Argument values as passed to __init__ will be saved, with
    the exception that if an instance attribute exists after __init__ has run that has
    the same name as an argument, the instance attribute will be saved instead of the
    argument value. This allows code within __init__ to process default arguments
    before they are saved.
    
    property_names is a dictionary {key:val}, where each val
        is a list [var1, var2, ...] of variables to be pulled from
        properties_dict and added to the property with name key (its location)
        
    internally they are all accessed by calling self.get_property()
    """
    def decorator(func):
        @wraps(func)
        def new_function(inst, *args, **kwargs):

            return_value = func(inst, *args, **kwargs)

            # Get a dict of the call arguments/keyword arguments by name:
            call_values = getcallargs(func, inst, *args, **kwargs)

            all_property_names = set()
            for names in property_names.values():
                all_property_names.update(names)

            property_values = {}
            for name in all_property_names:
                # If there is an instance attribute with that name, use that, otherwise
                # use the call value:
                if hasattr(inst, name):
                    property_values[name] = getattr(inst, name)
                else:
                    property_values[name] = call_values[name]

            # Save them:
            inst.set_properties(property_values, property_names)

            return return_value

        return new_function
    
    return decorator


class Device(object):
    description = 'Generic Device'
    allowed_children = None
    
    @set_passed_properties(
        property_names = {"device_properties": ["added_properties"]}
        )
    def __init__(self,name,parent_device,connection, call_parents_add_device=True, 
                 added_properties = {}, gui=None, worker=None, start_order=None, stop_order=None, **kwargs):
        # Verify that no invalid kwargs were passed and the set properties
        if len(kwargs) != 0:        
            raise LabscriptError('Invalid keyword arguments: %s.'%kwargs)

        if self.allowed_children is None:
            self.allowed_children = [Device]
        self.name = name
        self.parent_device = parent_device
        self.connection = connection
        self.start_order = start_order
        self.stop_order = stop_order
        if start_order is not None and not isinstance(start_order, int):
            raise TypeError("start_order must be int, not %s" % type(start_order).__name__)
        if stop_order is not None and not isinstance(stop_order, int):
            raise TypeError("stop_order must be int, not %s" % type(stop_order).__name__)
        self.child_devices = []
        
        # self._properties may be instantiated already
        if not hasattr(self, "_properties"):
            self._properties = {}
        for location in labscript_utils.properties.VALID_PROPERTY_LOCATIONS:
            if location not in self._properties:
                self._properties[location] = {}

        if parent_device and call_parents_add_device:
            # This is optional by keyword argument, so that subclasses
            # overriding __init__ can call call Device.__init__ early
            # on and only call self.parent_device.add_device(self)
            # a bit later, allowing for additional code in
            # between. If setting call_parents_add_device=False,
            # self.parent_device.add_device(self) *must* be called later
            # on, it is not optional.
            parent_device.add_device(self)
            
        # Check that the name doesn't already exist in the python namespace
        if name in locals() or name in globals() or name in _builtins_dict:
            raise LabscriptError('The device name %s already exists in the Python namespace. Please choose another.'%name)
        if name in keyword.kwlist:
            raise LabscriptError('%s is a reserved Python keyword.'%name +
                                 ' Please choose a different device name.')
                                     
        try:
            # Test that name is a valid Python variable name:
            exec('%s = None'%name)
            assert '.' not in name
        except:
            raise ValueError('%s is not a valid Python variable name.'%name)
        
        # Put self into the global namespace:
        _builtins_dict[name] = self
        
        # Add self to the compiler's device inventory
        compiler.inventory.append(self)
        
        # handle remote workers/gui interface
        if gui is not None or worker is not None:
            # remote GUI and worker
            if gui is not None:
                # if no worker is specified, assume it is the same as the gui
                if worker is None:
                    worker = gui
                    
                # check that worker and gui are appropriately typed
                if not isinstance(gui, _RemoteConnection):
                    raise LabscriptError('the "gui" argument for %s must be specified as a subclass of _RemoteConnection'%(self.name))
            else:
                # just remote worker
                gui = compiler._PrimaryBLACS
            
            if not isinstance(worker, _RemoteConnection):
                raise LabscriptError('the "worker" argument for %s must be specified as a subclass of _RemoteConnection'%(self.name))
            
            # check that worker is equal to, or a child of, gui
            if worker != gui and worker not in gui.get_all_children():
                print(gui.get_all_children())
                raise LabscriptError('The remote worker (%s) for %s must be a child of the specified gui (%s) '%(worker.name, self.name, gui.name))
                
            # store worker and gui as properties of the connection table
            self.set_property('gui', gui.name, 'connection_table_properties')
            self.set_property('worker', worker.name, 'connection_table_properties')
            
    
    # Method to set a property for this device. 
    #
    # The property will be stored in the connection table and be used 
    # during connection table comparisons.
    #
    # The value must satisfy eval(repr(value)) == value
    #
    # You cannot overwrite an existing property unless you set the 
    # overwrite flag to True on subsequent calls to this method
    # 
    # you can specify a location = "device_properties" or "connection_table_properties"
    # to set where these are stored.
    def set_property(self, name, value, location=None, overwrite=False):
        if location is None or location not in labscript_utils.properties.VALID_PROPERTY_LOCATIONS:
            raise LabscriptError('Device %s requests invalid property assignment %s for property %s'%(self.name, location, name))
            
        # if this try fails then self."location" may not be instantiated
        if not hasattr(self, "_properties"):
            self._properties = {}

        if location not in self._properties:
            self._properties[location] = {}

        selected_properties = self._properties[location]
        
        if name in selected_properties and not overwrite:
            raise LabscriptError('Device %s has had the property %s set more than once. This is not allowed unless the overwrite flag is explicitly set'%(self.name, name))

        selected_properties[name] = value

    def set_properties(self, properties_dict, property_names, overwrite = False):
        """
        Add one or a bunch of properties packed into properties_dict

        property_names is a dictionary {key:val, ...} where each val
            is a list [var1, var2, ...] of variables to be pulled from
            properties_dict and added to the property with name key (it's location)
        """
        for location, names in property_names.items():
            if not isinstance(names, list) and not isinstance(names, tuple):
                raise TypeError('%s names (%s) must be list or tuple, not %s'%(location, repr(names), str(type(names))))
            temp_dict = {key:val for key, val in properties_dict.items() if key in names}                  
            for (name, value) in temp_dict.items():
                self.set_property(name, value, 
                                  overwrite = overwrite, 
                                  location = location)
    
    # Method to get a property of this device already set using Device.set_property
    #
    # If the property is not already set, a default value will be returned 
    # if specified as the argument after 'name' if there is only one argument
    # after 'name' and the argument is either not a keyword argument or is a 
    # keyword argument with the name 'default'
    #
    # If the property is not already set, or one of the above conventions is
    # violated, a LabscriptError will be raised
    #
    # Example acceptable signatures:
    #
    # get_property('example')            # 'example will be returned if set, or an exception will be raised
    # get_property('example', 7)         # 7 will be returned if 'example' is not set
    # get_property('example', default=7) # 7 will be returned is 'example' is not set
    #
    #
    # Example signatures that WILL ALWAYS RAISE AN EXCEPTION:
    # get_property('example', 7, 8)
    # get_property('example', 7, default=9)
    # get_property('example', default=7, x=9)
    #
    # the named argument location, if passed, requests the keyword be searched
    # from only that location
    def get_property(self, name, location = None, *args, **kwargs):#default = None):
        if len(kwargs) == 1 and 'default' not in kwargs:
            raise LabscriptError('A call to %s.get_property had a keyword argument that was not name or default'%self.name)
        if len(args) + len(kwargs) > 1:
            raise LabscriptError('A call to %s.get_property has too many arguments and/or keyword arguments'%self.name)

        if (location is not None) and (location not in labscript_utils.properties.VALID_PROPERTY_LOCATIONS):
            raise LabscriptError('Device %s requests invalid property read location %s'%(self.name, location))
            
        # self._properties may not be instantiated
        if not hasattr(self, "_properties"):
            self._properties =  {}
        
        # Run through all keys of interest
        for key, val in self._properties.items():
            if (location is None or key == location) and (name in val):
               return val[name]
            
        if 'default' in kwargs:
            return kwargs['default']
        elif len(args) == 1:
            return args[0]
        else:
            raise LabscriptError('The property %s has not been set for device %s'%(name, self.name))

    def get_properties(self, location = None):
        """
        Get all properties in location
        
        If location is None we return all keys
        """
    
        # self._properties may not be instantiated
        if not hasattr(self, "_properties"):
            self._properties =  {}

        if location is not None:
            temp_dict = self._properties.get(location, {})
        else:
            temp_dict = {}
            for key,val in self._properties.items(): temp_dict.update(val)
                
        return temp_dict

    def add_device(self, device):
        if any([isinstance(device,DeviceClass) for DeviceClass in self.allowed_children]):
            self.child_devices.append(device)
        else:
            raise LabscriptError('Devices of type %s cannot be attached to devices of type %s.'%(device.description,self.description))
    
    @property    
    def pseudoclock_device(self):
        if isinstance(self, PseudoclockDevice):
            return self 
        parent = self.parent_device
        try:
            while parent is not None and not isinstance(parent,PseudoclockDevice):
                parent = parent.parent_device
            return parent
        except Exception as e:
            raise LabscriptError('Couldn\'t find parent pseudoclock device of %s, what\'s going on? Original error was %s.'%(self.name, str(e)))
    
    def quantise_to_pseudoclock(self, times):
        convert_back_to = None 
        if not isinstance(times, ndarray):
            if isinstance(times, list):
                convert_back_to = list
            elif isinstance(times, set):
                convert_back_to = set
            else:
                convert_back_to = float
            times = array(times)
        # quantise the times to the pseudoclock clock resolution
        times = (times/self.pseudoclock_device.clock_resolution).round()*self.pseudoclock_device.clock_resolution
        
        if convert_back_to is not None:
            times = convert_back_to(times)
        
        return times
    
    @property 
    def parent_clock_line(self):
        if isinstance(self, ClockLine):
            return self
        parent = self.parent_device
        try:
            while not isinstance(parent,ClockLine):
                parent = parent.parent_device
            return parent
        except Exception as e:
            raise LabscriptError('Couldn\'t find parent ClockLine of %s, what\'s going on? Original error was %s.'%(self.name, str(e)))
    
    @property
    def t0(self):
        """The earliest time output can be commanded from this device at the start of the experiment.
        This is nonzeo on secondary pseudoclock devices due to triggering delays."""
        parent = self.pseudoclock_device
        if parent is None or parent.is_master_pseudoclock:
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

    def init_device_group(self, hdf5_file):
        group = hdf5_file['/devices'].create_group(self.name)
        return group


class _PrimaryBLACS(Device):
    pass
    
class _RemoteConnection(Device):
    @set_passed_properties(
        property_names = {}
    )
    def __init__(self, name, parent=None, connection=None):
        if parent is None:
            # define a hidden parent of top level remote connections so that
            # "connection" is stored correctly in the connection table
            if compiler._PrimaryBLACS is None:
                compiler._PrimaryBLACS = _PrimaryBLACS('__PrimaryBLACS', None, None)
            parent = compiler._PrimaryBLACS
        Device.__init__(self, name, parent, connection)
        
        
class RemoteBLACS(_RemoteConnection):
    def __init__(self, name, host, port=7341, parent=None):
        _RemoteConnection.__init__(self, name, parent, "%s:%s"%(host, port))
        
        
class SecondaryControlSystem(_RemoteConnection):
    def __init__(self, name, host, port, parent=None):
        _RemoteConnection.__init__(self, name, parent, "%s:%s"%(host, port))
        

class IntermediateDevice(Device):
    
    @set_passed_properties(property_names = {})
    def __init__(self, name, parent_device, **kwargs):

        self.name = name
        # this should be checked here because it should only be connected a clockline
        # The allowed_children attribute of parent classes doesn't prevent this from being connected to something that accepts 
        # an instance of 'Device' as a child
        if parent_device is not None and not isinstance(parent_device, ClockLine):
            if not hasattr(parent_device, 'name'):
                parent_device_name = 'Unknown: not an instance of a labscript device class'
            else:
                parent_device_name = parent_device.name
            raise LabscriptError('Error instantiating device %s. The parent (%s) must be an instance of ClockLine.'%(name, parent_device_name))
        Device.__init__(self, name, parent_device, 'internal', **kwargs) # This 'internal' should perhaps be more descriptive?
 
  
class ClockLine(Device):
    description = 'Generic ClockLine'
    allowed_children = [IntermediateDevice]
    _clock_limit = None
    
    @set_passed_properties(property_names = {})
    def __init__(self, name, pseudoclock, connection, ramping_allowed = True, **kwargs):
        
        # TODO: Verify that connection is  valid connection of Pseudoclock.parent_device (the PseudoclockDevice)
        Device.__init__(self, name, pseudoclock, connection, **kwargs)
        self.ramping_allowed = ramping_allowed
        
    def add_device(self, device):
        Device.add_device(self, device)
        if getattr(device, 'clock_limit', None) is not None and (self._clock_limit is None or device.clock_limit < self.clock_limit):
            self._clock_limit = device.clock_limit
    
    # define a property to make sure no children overwrite this value themselves
    # The calculation of maximum clock_limit should be done by the add_device method above
    @property
    def clock_limit(self):
        # If no child device has specified a clock limit
        if self._clock_limit is None:
            # return the Pseudoclock clock_limit
            # TODO: Maybe raise an error instead?
            #       Maybe all Intermediate devices should be required to have a clock_limit?
            return self.parent_device.clock_limit
        return self._clock_limit

        
class Pseudoclock(Device):
    description = 'Generic Pseudoclock'
    allowed_children = [ClockLine]
    
    @set_passed_properties(property_names = {})
    def __init__(self, name, pseudoclock_device, connection, **kwargs):

        Device.__init__(self, name, pseudoclock_device, connection, **kwargs)
        self.clock_limit = pseudoclock_device.clock_limit
        self.clock_resolution = pseudoclock_device.clock_resolution
        
    def add_device(self, device):
        Device.add_device(self, device)
        #TODO: Maybe verify here that device.connection (the ClockLine connection) is a valid connection of the parent PseudoClockDevice
        #      Also see the same comment in ClockLine.__init__
        # if device.connection not in self.clock_lines:
            # self.clock_lines[

    def collect_change_times(self, all_outputs, outputs_by_clockline):
        """Asks all connected outputs for a list of times that they
        change state. Takes the union of all of these times. Note
        that at this point, a change from holding-a-constant-value
        to ramping-through-values is considered a single state
        change. The clocking times will be filled in later in the
        expand_change_times function, and the ramp values filled in with
        expand_timeseries."""
        change_times = {}
        all_change_times = []
        ramps_by_clockline = {}
        for clock_line, outputs in outputs_by_clockline.items():
            change_times.setdefault(clock_line, [])
            ramps_by_clockline.setdefault(clock_line, [])
            for output in outputs:
                # print 'output name: %s'%output.name
                output_change_times = output.get_change_times()
                # print output_change_times
                change_times[clock_line].extend(output_change_times)
                all_change_times.extend(output_change_times)
                ramps_by_clockline[clock_line].extend(output.get_ramp_times())
            
            # print 'initial_change_times for %s: %s'%(clock_line.name,change_times[clock_line])
        
        # Change to a set and back to get rid of duplicates:
        if not all_change_times:
            all_change_times.append(0)
        all_change_times.append(self.parent_device.stop_time)
        # include trigger times in change_times, so that pseudoclocks always have an instruction immediately following a wait:
        all_change_times.extend(self.parent_device.trigger_times)
        
        ####################################################################################################
        # Find out whether any other clockline has a change time during a ramp on another clockline.       #
        # If it does, we need to let the ramping clockline know it needs to break it's loop at that time   #
        ####################################################################################################
        # convert all_change_times to a numpy array
        all_change_times_numpy = array(all_change_times)
        
        # quantise the all change times to the pseudoclock clock resolution
        # all_change_times_numpy = (all_change_times_numpy/self.clock_resolution).round()*self.clock_resolution
        all_change_times_numpy = self.quantise_to_pseudoclock(all_change_times_numpy)
        
        # Loop through each clockline
        # print ramps_by_clockline
        for clock_line, ramps in ramps_by_clockline.items():
            # for each clockline, loop through the ramps on that clockline
            for ramp_start_time, ramp_end_time in ramps:
                # for each ramp, check to see if there is a change time in all_change_times which intersects
                # with the ramp. If there is, add a change time into this clockline at that point
                indices = np.where((ramp_start_time < all_change_times_numpy) & (all_change_times_numpy < ramp_end_time))
                for idx in indices[0]:
                    change_times[clock_line].append(all_change_times_numpy[idx])
                
        # Get rid of duplicates:
        all_change_times = list(set(all_change_times_numpy))
        all_change_times.sort()  
        
        # Check that the pseudoclock can handle updates this fast
        for i, t in enumerate(all_change_times[:-1]):
            dt = all_change_times[i+1] - t
            if dt < 1.0/self.clock_limit:
                raise LabscriptError('Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(all_change_times[i+1])) +
                                     'This Pseudoclock cannot support update delays shorter than %s sec.'%(str(1.0/self.clock_limit)))

        ####################################################################################################
        # For each clockline, make sure we have a change time for triggers, stop_time, t=0 and             #
        # check that no change tiems are too close together                                                #
        ####################################################################################################
        for clock_line, change_time_list in change_times.items():
            # include trigger times in change_times, so that pseudoclocks always have an instruction immediately following a wait:
            change_time_list.extend(self.parent_device.trigger_times)
            
            # If the device has no children, we still need it to have a
            # single instruction. So we'll add 0 as a change time:
            if not change_time_list:
                change_time_list.append(0)
            
            # quantise the all change times to the pseudoclock clock resolution
            # change_time_list = (array(change_time_list)/self.clock_resolution).round()*self.clock_resolution
            change_time_list = self.quantise_to_pseudoclock(change_time_list)
            
            # Get rid of duplicates if trigger times were already in the list:
            change_time_list = list(set(change_time_list))
            change_time_list.sort()

            # index to keep track of in all_change_times
            j = 0
            # Check that no two instructions are too close together:
            for i, t in enumerate(change_time_list[:-1]):
                dt = change_time_list[i+1] - t
                if dt < 1.0/clock_line.clock_limit:
                    raise LabscriptError('Commands have been issued to devices attached to clockline %s at t= %s s and %s s. '%(clock_line.name, str(t),str(change_time_list[i+1])) +
                                         'One or more connected devices on ClockLine %s cannot support update delays shorter than %s sec.'%(clock_line.name, str(1.0/clock_line.clock_limit)))
                                         
                # increment j until we reach the current time
                while all_change_times[j] < t and j < len(all_change_times)-1:
                    j += 1
                # j should now index all_change_times at "t"
                # Check that the next all change_time is not too close (and thus would force this clock tick to be faster than the clock_limit)
                dt = all_change_times[j+1] - t
                if dt < 1.0/clock_line.clock_limit:
                    raise LabscriptError('Commands have been issued to devices attached to %s at t= %s s and %s s. '%(self.name, str(t),str(all_change_times[j+1])) +
                                         'One or more connected devices on ClockLine %s cannot support update delays shorter than %s sec.'%(clock_line.name, str(1.0/clock_line.clock_limit)))
            
            # Also add the stop time as as change time. First check that it isn't too close to the time of the last instruction:
            if not self.parent_device.stop_time in change_time_list:
                dt = self.parent_device.stop_time - change_time_list[-1]
                if abs(dt) < 1.0/clock_line.clock_limit:
                    raise LabscriptError('The stop time of the experiment is t= %s s, but the last instruction for a device attached to %s is at t= %s s. '%( str(self.stop_time), self.name, str(change_time_list[-1])) +
                                         'One or more connected devices cannot support update delays shorter than %s sec. Please set the stop_time a bit later.'%str(1.0/clock_line.clock_limit))
                
                change_time_list.append(self.parent_device.stop_time)

            # Sort change times so self.stop_time will be in the middle
            # somewhere if it is prior to the last actual instruction. Whilst
            # this means the user has set stop_time in error, not catching
            # the error here allows it to be caught later by the specific
            # device that has more instructions after self.stop_time. Thus
            # we provide the user with sligtly more detailed error info.
            change_time_list.sort()
            
            # because we made the list into a set and back to a list, it is now a different object
            # so modifying it won't update the list in the dictionary.
            # So store the updated list in the dictionary
            change_times[clock_line] = change_time_list
        return all_change_times, change_times
    
    def expand_change_times(self, all_change_times, change_times, outputs_by_clockline):
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
        # for output in outputs:
            # if output.parent_device.clock_type != 'slow clock':            
                # if output.parent_device.clock_type not in all_times:
                    # all_times[output.parent_device.clock_type] = []
                # if output.parent_device.clock_type not in clocks_in_use:
                    # clocks_in_use.append(output.parent_device.clock_type)
        
        clock = []
        clock_line_current_indices = {}
        for clock_line, outputs in outputs_by_clockline.items():
            clock_line_current_indices[clock_line] = 0
            all_times[clock_line] = []
        
        # iterate over all change times
        # for clock_line, times in change_times.items():
            # print '%s: %s'%(clock_line.name, times)
        for i, time in enumerate(all_change_times):
            if time in self.parent_device.trigger_times[1:]:
                # A wait instruction:
                clock.append('WAIT')
                
            # list of enabled clocks
            enabled_clocks = []
            enabled_looping_clocks = []
            # enabled_non_looping_clocks = []
            
            # update clock_line indices
            for clock_line in clock_line_current_indices:
                try:
                    while change_times[clock_line][clock_line_current_indices[clock_line]] < time:
                        clock_line_current_indices[clock_line] += 1
                except IndexError:
                    # Fix the index to the last one
                    clock_line_current_indices[clock_line] = len(change_times[clock_line]) - 1
                    # print a warning
                    message = ''.join(['WARNING: ClockLine %s has it\'s last change time at t=%.15f but another ClockLine has a change time at t=%.15f. '%(clock_line.name, change_times[clock_line][-1], time), 
                              'This should never happen, as the last change time should always be the time passed to stop(). ', 
                              'Perhaps you have an instruction after the stop time of the experiment?'])
                    sys.stderr.write(message+'\n')
                    
                # Let's work out which clock_lines are enabled for this instruction
                if time == change_times[clock_line][clock_line_current_indices[clock_line]]:
                    enabled_clocks.append(clock_line)
            
            # what's the fastest clock rate?
            maxrate = 0
            local_clock_limit = self.clock_limit # the Pseudoclock clock limit
            for clock_line in enabled_clocks:
                for output in outputs_by_clockline[clock_line]:
                    # Check if output is sweeping and has highest clock rate
                    # so far. If so, store its clock rate to max_rate:
                    if hasattr(output,'timeseries') and isinstance(output.timeseries[clock_line_current_indices[clock_line]],dict):
                        if clock_line not in enabled_looping_clocks:
                            enabled_looping_clocks.append(clock_line)
                                
                        if output.timeseries[clock_line_current_indices[clock_line]]['clock rate'] > maxrate:
                            # It does have the highest clock rate? Then store that rate to max_rate:
                            maxrate = output.timeseries[clock_line_current_indices[clock_line]]['clock rate']
                
                        # only check this for ramping clock_lines
                        # non-ramping clock-lines have already had the clock_limit checked within collect_change_times()
                        if local_clock_limit > clock_line.clock_limit:
                            local_clock_limit = clock_line.clock_limit
                        
            # find non-looping clocks
            # for clock_line in enabled_clocks:
                # if clock_line not in enabled_looping_clocks:
                    # enabled_non_looping_clocks.append(clock_line)
            
            if maxrate:
                # round to the nearest clock rate that the pseudoclock can actually support:
                period = 1/maxrate
                quantised_period = period/self.clock_resolution
                quantised_period = round(quantised_period)
                period = quantised_period*self.clock_resolution
                maxrate = 1/period
            if maxrate > local_clock_limit:
                raise LabscriptError('At t = %s sec, a clock rate of %s Hz was requested. '%(str(time),str(maxrate)) + 
                                    'One or more devices connected to %s cannot support clock rates higher than %sHz.'%(str(self.name),str(local_clock_limit)))
                
            if maxrate:
                # If there was ramping at this timestep, how many clock ticks fit before the next instruction?
                n_ticks, remainder = divmod((all_change_times[i+1] - time)*maxrate,1)
                n_ticks = int(n_ticks)
                # Can we squeeze the final clock cycle in at the end?
                if remainder and remainder/float(maxrate) >= 1/float(local_clock_limit):
                    # Yes we can. Clock speed will be as
                    # requested. Otherwise the final clock cycle will
                    # be too long, by the fraction 'remainder'.
                    n_ticks += 1
                duration = n_ticks/float(maxrate) # avoiding integer division
                ticks = linspace(time,time + duration,n_ticks,endpoint=False)
                
                for clock_line in enabled_clocks:
                    if clock_line in enabled_looping_clocks:
                        all_times[clock_line].append(ticks)
                    else:
                        all_times[clock_line].append(time)
                
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
                        clock.append({'start': time, 'reps': 1, 'step': 1/float(maxrate), 'enabled_clocks':enabled_clocks})
                        clock.append({'start': time + 1/float(maxrate), 'reps': n_ticks-2, 'step': 1/float(maxrate), 'enabled_clocks':enabled_looping_clocks})
                    else:
                        clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate), 'enabled_clocks':enabled_clocks})
                        
                    # clock.append({'start': time, 'reps': n_ticks-1, 'step': 1/float(maxrate), 'enabled_clocks':enabled_clocks})
                # The last clock tick has a different duration depending on the next step. 
                clock.append({'start': ticks[-1], 'reps': 1, 'step': all_change_times[i+1] - ticks[-1], 'enabled_clocks':enabled_clocks if n_ticks == 1 else enabled_looping_clocks})
            else:
                for clock_line in enabled_clocks:
                    all_times[clock_line].append(time)
                    
                try: 
                    # If there was no ramping, here is a single clock tick:
                    clock.append({'start': time, 'reps': 1, 'step': all_change_times[i+1] - time, 'enabled_clocks':enabled_clocks})
                except IndexError:
                    if i != len(all_change_times) - 1:
                        raise
                    if self.parent_device.stop_time > time:
                        # There is no next instruction. Hold the last clock
                        # tick until self.parent_device.stop_time.
                        raise Exception('This shouldn\'t happen -- stop_time should always be equal to the time of the last instruction. Please report a bug.')
                        # I commented this out because it is after a raised exception so never ran.... - Phil
                        # clock.append({'start': time, 'reps': 1, 'step': self.parent_device.stop_time - time,'slow_clock_tick':True}) 
                    # Error if self.parent_device.stop_time has been set to less
                    # than the time of the last instruction:
                    elif self.parent_device.stop_time < time:
                        raise LabscriptError('%s %s has more instructions (at t=%.15f) after the experiment\'s stop time (t=%.15f).'%(self.description,self.name, time, self.parent_device.stop_time))
                    # If self.parent_device.stop_time is the same as the time of the last
                    # instruction, then we'll get the last instruction
                    # out still, so that the total number of clock
                    # ticks matches the number of data points in the
                    # Output.raw_output arrays. We'll make this last
                    # cycle be at ten times the minimum step duration.
                    else:
                        # find the slowest clock_limit
                        enabled_clocks = []
                        local_clock_limit = 1.0/self.clock_resolution # the Pseudoclock clock limit
                        for clock_line, outputs in outputs_by_clockline.items():
                            if local_clock_limit > clock_line.clock_limit:
                                local_clock_limit = clock_line.clock_limit
                            enabled_clocks.append(clock_line)
                        clock.append({'start': time, 'reps': 1, 'step': 10.0/self.clock_limit, 'enabled_clocks':enabled_clocks})
        # for row in clock:
            # print row
        return all_times, clock
    
    def get_outputs_by_clockline(self):
        outputs_by_clockline = {}
        for clock_line in self.child_devices:
            if isinstance(clock_line, ClockLine):
                outputs_by_clockline[clock_line] = []

        all_outputs = self.get_all_outputs()
        for output in all_outputs:
            # TODO: Make this a bit more robust (can we assume things always have this hierarchy?)
            clock_line = output.parent_clock_line
            assert clock_line.parent_device == self
            outputs_by_clockline[clock_line].append(output)

        return all_outputs, outputs_by_clockline

    def generate_clock(self):
        all_outputs, outputs_by_clockline = self.get_outputs_by_clockline()
        
        # Get change_times for all outputs, and also grouped by clockline
        all_change_times, change_times = self.collect_change_times(all_outputs, outputs_by_clockline)
               
        # for each clock line
        for clock_line, clock_line_change_times in change_times.items():
            # and for each output on the clockline
            for output in outputs_by_clockline[clock_line]:
                # call make_timeseries to expand the list of instructions for each change_time on this clock line
                output.make_timeseries(clock_line_change_times)

        # now generate the clock meta data for the Pseudoclock
        # also generate everytime point each clock line will tick (expand ramps)
        all_times, self.clock = self.expand_change_times(all_change_times, change_times, outputs_by_clockline)
        
        # for each clockline
        for clock_line, outputs in outputs_by_clockline.items():
            # and for each output
            for output in outputs:
                # evaluate the output at each time point the clock line will tick at
                output.expand_timeseries(all_times[clock_line])
                
        # TODO: is this needed? Let's say no...
        # self.all_change_times = fastflatten(all_change_times, float)
        
        # Flatten the clock line times for use by the child devices for writing instruction tables
        # TODO: (if this needed or was it just for runviewer meta data that we don't need anymore?)
        self.times = {}
        for clock_line, time_array in all_times.items():
            self.times[clock_line] = fastflatten(time_array,float)
        
    def generate_code(self, hdf5_file):
        self.generate_clock()
        Device.generate_code(self, hdf5_file)
        

class TriggerableDevice(Device):
    trigger_edge_type = 'rising'
    minimum_recovery_time = 0
    # A class devices should inherit if they do
    # not require a pseudoclock, but do require a trigger.
    # This enables them to have a Trigger device as a parent
    
    @set_passed_properties(property_names = {})
    def __init__(self, name, parent_device, connection, parentless=False, **kwargs):

        if None in [parent_device, connection] and not parentless:
            raise LabscriptError('No parent specified. If this device does not require a parent, set parentless=True')
        if isinstance(parent_device, Trigger):
            if self.trigger_edge_type != parent_device.trigger_edge_type:
                raise LabscriptError('Trigger edge type for %s is \'%s\', ' % (name, self.trigger_edge_type) + 
                                      'but existing Trigger object %s ' % parent_device.name +
                                      'has edge type \'%s\'' % parent_device.trigger_edge_type)
            self.trigger_device = parent_device
        elif parent_device is not None:
            # Instantiate a trigger object to be our parent:
            self.trigger_device = Trigger(name + '_trigger', parent_device, connection, self.trigger_edge_type)
            parent_device = self.trigger_device
            connection = 'trigger'
            
        self.__triggers = []
        Device.__init__(self, name, parent_device, connection, **kwargs)

    def trigger(self, t, duration):
        """Request parent trigger device to produce a trigger at time t with given
        duration."""
        # Only ask for a trigger if one has not already been requested by another device
        # attached to the same trigger:
        already_requested = False
        for other_device in self.trigger_device.child_devices:
            if other_device is not self:
                for other_t, other_duration in other_device.__triggers:
                    if t == other_t and duration == other_duration:
                        already_requested = True
        if not already_requested:
            self.trigger_device.trigger(t, duration)

        # Check for triggers too close together (check for overlapping triggers already
        # performed in Trigger.trigger()):
        start = t
        end = t + duration
        for other_t, other_duration in self.__triggers:
            other_start = other_t
            other_end = other_t + other_duration
            if (
                abs(other_start - end) < self.minimum_recovery_time
                or abs(other_end - start) < self.minimum_recovery_time
            ):
                msg = """%s %s has two triggers closer together than the minimum
                    recovery time: one at t = %fs for %fs, and another at t = %fs for
                    %fs. The minimum recovery time is %fs."""
                msg = msg % (
                    self.description,
                    self.name,
                    t,
                    duration,
                    start,
                    duration,
                    self.minimum_recovery_time,
                )
                raise ValueError(dedent(msg))

        self.__triggers.append([t, duration])

    def do_checks(self):
        # Check that all devices sharing a trigger device have triggers when we have triggers:
        for device in self.trigger_device.child_devices:
            if device is not self:
                for trigger in self.__triggers:
                    if trigger not in device.__triggers:
                        start, duration = trigger
                        raise LabscriptError('TriggerableDevices %s and %s share a trigger. ' % (self.name, device.name) + 
                                             '%s has a trigger at %fs for %fs, ' % (self.name, start, duration) +
                                             'but there is no matching trigger for %s. ' % device.name +
                                             'Devices sharing a trigger must have identical trigger times and durations.')

    def generate_code(self, hdf5_file):
        self.do_checks()
        Device.generate_code(self, hdf5_file)

    
class PseudoclockDevice(TriggerableDevice):
    description = 'Generic Pseudoclock Device'
    allowed_children = [Pseudoclock]
    trigger_edge_type = 'rising'
    # How long after a trigger the next instruction is actually output:
    trigger_delay = 0
    # How long a trigger line must remain high/low in order to be detected:
    trigger_minimum_duration = 0 
    # How long after the start of a wait instruction the device is actually capable of resuming:
    wait_delay = 0
    
    @set_passed_properties(property_names = {})
    def __init__(self, name, trigger_device=None, trigger_connection=None, **kwargs):

        if trigger_device is None:
            for device in compiler.inventory:
                if isinstance(device, PseudoclockDevice) and device.is_master_pseudoclock:
                    raise LabscriptError('There is already a master pseudoclock device: %s.'%device.name + 
                                         'There cannot be multiple master pseudoclock devices - please provide a trigger_device for one of them.')
            TriggerableDevice.__init__(self, name, parent_device=None, connection=None, parentless=True, **kwargs)
        else:
            # The parent device declared was a digital output channel: the following will
            # automatically instantiate a Trigger for it and set it as self.trigger_device:
            TriggerableDevice.__init__(self, name, parent_device=trigger_device, connection=trigger_connection, **kwargs)
            # Ensure that the parent pseudoclock device is, in fact, the master pseudoclock device.
            if not self.trigger_device.pseudoclock_device.is_master_pseudoclock:
                raise LabscriptError('All secondary pseudoclock devices must be triggered by a device being clocked by the master pseudoclock device.' +
                                     'Pseudoclocks triggering each other in series is not supported.')
        self.trigger_times = []
        self.wait_times = []
        self.initial_trigger_time = 0
    
    @property    
    def is_master_pseudoclock(self):
        return self.parent_device is None
    
    def set_initial_trigger_time(self, t):
        t = round(t,10)
        if compiler.start_called:
            raise LabscriptError('Initial trigger times must be set prior to calling start()')
        if self.is_master_pseudoclock:
            raise LabscriptError('Initial trigger time of master clock is always zero, it cannot be changed.')
        else:
            self.initial_trigger_time = t
            
    def trigger(self, t, duration, wait_delay = 0):
        """Ask the trigger device to produce a digital pulse of a given duration to trigger this pseudoclock"""
        if type(t) in [str, bytes] and t == 'initial':
            t = self.initial_trigger_time
        t = round(t,10)
        if self.is_master_pseudoclock:
            if compiler.wait_monitor is not None:
                # Make the wait monitor pulse to signify starting or resumption of the experiment:
                compiler.wait_monitor.trigger(t, duration)
            self.trigger_times.append(t)
        else:
            TriggerableDevice.trigger(self, t, duration)
            self.trigger_times.append(round(t + wait_delay,10))
            
    def do_checks(self, outputs):
        """Basic error checking to ensure the user's instructions make sense"""
        for output in outputs:
            output.do_checks(self.trigger_times)
            
    def offset_instructions_from_trigger(self, outputs):
        for output in outputs:
            output.offset_instructions_from_trigger(self.trigger_times)
        
        if not self.is_master_pseudoclock:
            # Store the unmodified initial_trigger_time
            initial_trigger_time = self.trigger_times[0]
            # Adjust the stop time relative to the last trigger time
            self.stop_time = round(self.stop_time - initial_trigger_time - self.trigger_delay * len(self.trigger_times),10)
            
            # Modify the trigger times themselves so that we insert wait instructions at the right times:
            self.trigger_times = [round(t - initial_trigger_time - i*self.trigger_delay,10) for i, t in enumerate(self.trigger_times)]
        
        # quantise the trigger times and stop time to the pseudoclock clock resolution
        self.trigger_times = self.quantise_to_pseudoclock(self.trigger_times)
        self.stop_time = self.quantise_to_pseudoclock(self.stop_time)
                            
    def generate_code(self, hdf5_file):
        outputs = self.get_all_outputs()
        self.do_checks(outputs)
        self.offset_instructions_from_trigger(outputs)
        Device.generate_code(self, hdf5_file)
        
    
class Output(Device):
    description = 'generic output'
    allowed_states = {}
    dtype = float64
    scale_factor = 1
    
    @set_passed_properties(property_names = {})
    def __init__(self,name,parent_device,connection,limits = None,unit_conversion_class = None, unit_conversion_parameters = None, default_value=None, **kwargs):
        Device.__init__(self,name,parent_device,connection, **kwargs)

        self.instructions = {}
        self.ramp_limits = [] # For checking ramps don't overlap
        if default_value is not None:
            self.default_value = default_value
        if not unit_conversion_parameters:
            unit_conversion_parameters = {}
        self.unit_conversion_class = unit_conversion_class
        self.set_properties(unit_conversion_parameters,
                            {'unit_conversion_parameters': list(unit_conversion_parameters.keys())})

        # Instantiate the calibration
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
            if not isinstance(limits,tuple) or len(limits) != 2:
                raise LabscriptError('The limits for "%s" must be tuple of length 2. Eg. limits=(1,2)'%(self.name))
            if limits[0] > limits[1]:
                raise LabscriptError('The first element of the tuple must be lower than the second element. Eg limits=(1,2), NOT limits=(2,1)')
        # Save limits even if they are None        
        self.limits = limits
    
    @property
    def clock_limit(self):
        parent = self.parent_clock_line
        return parent.clock_limit
    
    @property
    def trigger_delay(self):
        """The earliest time output can be commanded from this device after a trigger.
        This is nonzeo on secondary pseudoclocks due to triggering delays."""
        parent = self.pseudoclock_device
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
        delay = compiler.wait_delay if self.pseudoclock_device.is_master_pseudoclock else 0
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
                 'Due to the delay in triggering its pseudoclock device, the earliest output possible is at t=%s.'%str(self.t0)])
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
            if not self.parent_clock_line.ramping_allowed:
                raise LabscriptError('%s %s is on clockline that does not support ramping. '%(self.description, self.name) + 
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
                if (instruction < self.limits[0]) or (instruction > self.limits[1]):
                    raise LabscriptError('You cannot program the value %s (base units) to %s as it falls outside the limits (%d to %d)'%(str(instruction), self.name, self.limits[0], self.limits[1]))
        self.instructions[time] = instruction
    
    def do_checks(self, trigger_times):
        """Basic error checking to ensure the user's instructions make sense"""
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
        for instruction in list(self.instructions.values()):
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
            offset = round(offset,10)
            if isinstance(instruction,dict):
                offset_instruction = instruction.copy()
                offset_instruction['end time'] = self.quantise_to_pseudoclock(round(instruction['end time'] - offset,10))
                offset_instruction['initial time'] = self.quantise_to_pseudoclock(round(instruction['initial time'] - offset,10))
            else:
                offset_instruction = instruction
                
            offset_instructions[self.quantise_to_pseudoclock(round(t - offset,10))] = offset_instruction
        self.instructions = offset_instructions
            
        # offset each of the ramp_limits for use in the calculation within Pseudoclock/ClockLine
        # so that the times in list are consistent with the ones in self.instructions
        for i, times in enumerate(self.ramp_limits):
            n_triggers_prior = len([time for time in trigger_times if time < times[0]])
            # The cumulative offset at this point in time:
            offset = self.trigger_delay * n_triggers_prior + trigger_times[0]
            offset = round(offset,10)
            
            # offset start and end time of ramps
            # NOTE: This assumes ramps cannot proceed across a trigger command
            #       (for instance you cannot ramp an output across a WAIT)
            self.ramp_limits[i] = (self.quantise_to_pseudoclock(round(times[0]-offset,10)), self.quantise_to_pseudoclock(round(times[1]-offset,10)))
            
    def get_change_times(self):
        """If this function is being called, it means that the parent
        Pseudoclock has requested a list of times that this output changes
        state."""        
        times = list(self.instructions.keys())
        times.sort()

        current_dict_time = None
        for time in times:
            if isinstance(self.instructions[time], dict) and current_dict_time is None:
                current_dict_time = self.instructions[time]
            elif current_dict_time is not None and current_dict_time['initial time'] < time < current_dict_time['end time']:
                err = ("{:s} {:s} has an instruction at t={:.10f}s. This instruction collides with a ramp on this output at that time. ".format(self.description, self.name, time)+
                       "The collision {:s} is happening from {:.10f}s till {:.10f}s".format(current_dict_time['description'], current_dict_time['initial time'], current_dict_time['end time']))
                raise LabscriptError(err)

        self.times = times
        return times
        
    def get_ramp_times(self):
        return self.ramp_limits
    
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
        # If this output is not ramping, then its timeseries should
        # not be expanded. It's already as expanded as it'll get.
        if not self.parent_clock_line.ramping_allowed:
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
                            raise LabscriptError('The function %s called on "%s" at t=%d generated a value which falls outside the base unit limits (%d to %d)'%(self.timeseries[i]['function'],self.name,midpoints[0],self.limits[0],self.limits[1]))
                else:
                    outarray = empty(len(time),dtype=self.dtype)
                    outarray.fill(self.timeseries[i])
                outputarray.append(outarray)
            else:
                outputarray.append(self.timeseries[i])
        del self.timeseries # don't need this any more.
        self.raw_output = fastflatten(outputarray, self.dtype)
        

class AnalogQuantity(Output):
    description = 'analog quantity'
    default_value = 0

    def _check_truncation(self, truncation, min=0, max=1):
        if not (min <= truncation <= max):
            raise LabscriptError(
                'Truncation argument must be between %f and %f (inclusive), but is %f.' % (min, max, truncation))

    def ramp(self, t, duration, initial, final, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            # if start and end value are the same, we don't need to ramp and can save the sample ticks etc
            if initial == final:
                self.constant(t, initial, units)
                if not config.suppress_mild_warnings and not config.suppress_all_warnings:
                    message = ''.join(['WARNING: AnalogOutput \'%s\' has the same initial and final value at time t=%.10fs with duration %.10fs. In order to save samples and clock ticks this instruction is replaced with a constant output. '%(self.name, t, duration)])
                    sys.stderr.write(message + '\n')
            else:
                self.add_instruction(t, {'function': functions.ramp(round(t + duration, 10) - round(t, 10), initial, final), 'description': 'linear ramp',
                                     'initial time': t, 'end time': t + truncation * duration, 'clock rate': samplerate, 'units': units})
        return truncation * duration

    def sine(self, t, duration, amplitude, angfreq, phase, dc_offset, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            self.add_instruction(t, {'function': functions.sine(round(t + duration, 10) - round(t, 10), amplitude, angfreq, phase, dc_offset), 'description': 'sine wave',
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def sine_ramp(self, t, duration, initial, final, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            self.add_instruction(t, {'function': functions.sine_ramp(round(t + duration, 10) - round(t, 10), initial, final), 'description': 'sinusoidal ramp',
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def sine4_ramp(self, t, duration, initial, final, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            self.add_instruction(t, {'function': functions.sine4_ramp(round(t + duration, 10) - round(t, 10), initial, final), 'description': 'sinusoidal ramp',
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def sine4_reverse_ramp(self, t, duration, initial, final, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            self.add_instruction(t, {'function': functions.sine4_reverse_ramp(round(t + duration, 10) - round(t, 10), initial, final), 'description': 'sinusoidal ramp',
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def exp_ramp(self, t, duration, initial, final, samplerate, zero=0, units=None, truncation=None, truncation_type='linear', **kwargs):
        """Exponential ramp whose rate of change is set by an asymptotic value (zero argument).
        
        Parameters:
        t, duration : time to start the ramp and its duration
        initial, final : initial and final values of the ramp (sans truncation)
        zero: asymptotic value of the exponential decay/rise, i.e. limit as t --> inf
        samplerate: rate to sample the function
        units: unit conversion to apply to specified values before generating raw output
        truncation_type: 'linear' or 'exponential'
            'linear' truncation stops the ramp when it reaches the value given by the 
            truncation parameter, which must be between initial and final
            'exponential' truncation stops the ramp after a period of truncation*duration
            In this instance, the truncation parameter should be between 0 (full truncation)
            and 1 (no truncation).   
        """
        # Backwards compatibility for old kwarg names
        if 'trunc' in kwargs:
            truncation = kwargs.pop('trunc')
        if 'trunc_type' in kwargs:
            truncation_type = kwargs.pop('trunc_type')
        if truncation is not None:
            # Computed the truncated duration based on the truncation_type
            if truncation_type == 'linear':
                self._check_truncation(truncation, min(initial, final), max(initial, final))
                # Truncate the ramp when it reaches the value truncation
                trunc_duration = duration * \
                    log((initial-zero)/(truncation-zero)) / \
                    log((initial-zero)/(final-zero))
            elif truncation_type == 'exponential':
                # Truncate the ramps duration by a fraction truncation
                self._check_truncation(truncation)
                trunc_duration = truncation * duration
            else:
                raise LabscriptError(
                    'Truncation type for exp_ramp not supported. Must be either linear or exponential.')
        else:
            trunc_duration = duration
        if trunc_duration > 0:
            self.add_instruction(t, {'function': functions.exp_ramp(round(t + duration, 10) - round(t, 10), initial, final, zero), 'description': 'exponential ramp',
                                     'initial time': t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
        return trunc_duration

    def exp_ramp_t(self, t, duration, initial, final, time_constant, samplerate, units=None, truncation=None, truncation_type='linear', **kwargs):
        """Exponential ramp whose rate of change is set by the time_constant.
        
        Parameters:
        t, duration : time to start the ramp and its duration
        initial, final : initial and final values of the ramp (sans truncation)
        time_constant: 1/e time of the exponential decay/rise
        samplerate: rate to sample the function
        units: unit conversion to apply to specified values before generating raw output
        truncation_type: 'linear' or 'exponential'
            'linear' truncation stops the ramp when it reaches the value given by the 
            truncation parameter, which must be between initial and final
            'exponential' truncation stops the ramp after a period of truncation*duration
            In this instance, the truncation parameter should be between 0 (full truncation)
            and 1 (no truncation).   
        """
        # Backwards compatibility for old kwarg names
        if 'trunc' in kwargs:
            truncation = kwargs.pop('trunc')
        if 'trunc_type' in kwargs:
            truncation_type = kwargs.pop('trunc_type')
        if truncation is not None:
            zero = (final-initial*exp(-duration/time_constant)) / \
                (1-exp(-duration/time_constant))
            if truncation_type == 'linear':
                self._check_truncation(truncation, min(initial, final), max(initial, final))
                trunc_duration = time_constant * \
                    log((initial-zero)/(truncation-zero))
            elif truncation_type == 'exponential':
                self._check_truncation(truncation)
                trunc_duration = truncation * duration
            else:
                raise LabscriptError(
                    'Truncation type for exp_ramp_t not supported. Must be either linear or exponential.')
        else:
            trunc_duration = duration
        if trunc_duration > 0:
            self.add_instruction(t, {'function': functions.exp_ramp_t(round(t + duration, 10) - round(t, 10), initial, final, time_constant), 'description': 'exponential ramp with time consntant',
                                     'initial time': t, 'end time': t + trunc_duration, 'clock rate': samplerate, 'units': units})
        return trunc_duration

    def piecewise_accel_ramp(self, t, duration, initial, final, samplerate, units=None, truncation=1.):
        self._check_truncation(truncation)
        if truncation > 0:
            self.add_instruction(t, {'function': functions.piecewise_accel(round(t + duration, 10) - round(t, 10), initial, final), 'description': 'piecewise linear accelleration ramp',
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def customramp(self, t, duration, function, *args, **kwargs):
        units = kwargs.pop('units', None)
        samplerate = kwargs.pop('samplerate')
        truncation = kwargs.pop('truncation', 1.)
        self._check_truncation(truncation)

        def custom_ramp_func(t_rel):
            """The function that will return the result of the user's function,
            evaluated at relative times t_rel from 0 to duration"""
            return function(t_rel, round(t + duration, 10) - round(t, 10), *args, **kwargs)

        if truncation > 0:
            self.add_instruction(t, {'function': custom_ramp_func, 'description': 'custom ramp: %s' % function.__name__,
                                     'initial time': t, 'end time': t + truncation*duration, 'clock rate': samplerate, 'units': units})
        return truncation*duration

    def constant(self,t,value,units=None):
        # verify that value can be converted to float
        try:
            val = float(value)
        except:
            raise LabscriptError('in constant, value cannot be converted to float')
        self.add_instruction(t, value, units)
        
      
class AnalogOut(AnalogQuantity):
    description = 'analog output'
    
    
class StaticAnalogQuantity(Output):
    description = 'static analog quantity'
    default_value = 0.0
    
    @set_passed_properties(property_names = {})
    def __init__(self, *args, **kwargs):
        Output.__init__(self, *args, **kwargs)
        self._static_value = None
        
    def constant(self, value, units=None):
        if self._static_value is None:
            # If we have units specified, convert the value
            if units is not None:
                # Apply the unit calibration now
                value = self.apply_calibration(value, units)
            # if we have limits, check the value is valid
            if self.limits:
                minval, maxval = self.limits
                if not minval <= value <= maxval:
                    raise LabscriptError('You cannot program the value %s (base units) to %s as it falls outside the limits (%s to %s)'%(str(value), self.name, str(self.limits[0]), str(self.limits[1])))
            self._static_value = value
        else:
            raise LabscriptError('%s %s has already been set to %s (base units). It cannot also be set to %s (%s).'%(self.description, self.name, str(self._static_value), str(value),units if units is not None else "base units"))
    
    def get_change_times(self):
        return []  # Return an empty list as the calling function at the pseudoclock level expects a list
        
    def make_timeseries(self,change_times):
        pass
    
    def expand_timeseries(self,*args,**kwargs):
        self.raw_output = array([self.static_value], dtype=self.dtype)

    @property
    def static_value(self):
        if self._static_value is None:
            if not config.suppress_mild_warnings and not config.suppress_all_warnings:
                sys.stderr.write(' '.join(['WARNING:', self.name, 'has no value set. It will be set to %s.\n'%self.instruction_to_string(self.default_value)]))
            self._static_value = self.default_value
        return self._static_value
        
class StaticAnalogOut(StaticAnalogQuantity):
    description = 'static analog output'
        
class DigitalQuantity(Output):
    description = 'digital quantity'
    allowed_states = {1:'high', 0:'low'}
    default_value = 0
    dtype = uint32
    
    # Redefine __init__ so that you cannot define a limit or calibration for DO
    @set_passed_properties(property_names = {"connection_table_properties": ["inverted"]})
    def __init__(self, name, parent_device, connection, inverted=False, **kwargs):
        Output.__init__(self,name,parent_device,connection, **kwargs)
        self.inverted = bool(inverted)

    def go_high(self,t):
        self.add_instruction(t, 1)

    def go_low(self,t):
        self.add_instruction(t, 0)

    def enable(self,t):
        if self.inverted:
            self.go_low(t)
        else:
            self.go_high(t)

    def disable(self,t):
        if self.inverted:
            self.go_high(t)
        else:
            self.go_low(t)

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
    default_value = 0
    
    @set_passed_properties(property_names = {})
    def __init__(self, *args, **kwargs):
        DigitalQuantity.__init__(self, *args, **kwargs)
        self._static_value = None
        
    def go_high(self):
        if self._static_value is None:
            self.add_instruction(0,1)
            self._static_value = 1
        else:
            raise LabscriptError('%s %s has already been set to %s. It cannot also be set to %s.'%(self.description, self.name, self.instruction_to_string[self._static_value], self.instruction_to_string[value]))
            
    def go_low(self):
        if self._static_value is None:
            self.add_instruction(0,0) 
            self._static_value = 0
        else:
            raise LabscriptError('%s %s has already been set to %s. It cannot also be set to %s.'%(self.description, self.name, self.instruction_to_string[self._static_value], self.instruction_to_string[value]))
    
    def get_change_times(self):
        return []  # Return an empty list as the calling function at the pseudoclock level expects a list
    
    def make_timeseries(self,change_times):
        pass
    
    def expand_timeseries(self,*args,**kwargs):
        self.raw_output = array([self.static_value], dtype=self.dtype)

    @property
    def static_value(self):
        if self._static_value is None:
            if not config.suppress_mild_warnings and not config.suppress_all_warnings:
                sys.stderr.write(' '.join(['WARNING:', self.name, 'has no value set. It will be set to %s.\n'%self.instruction_to_string(self.default_value)]))
            self._static_value = self.default_value
        return self._static_value
    

class StaticDigitalOut(StaticDigitalQuantity):
    description = 'static digital output'
        
class AnalogIn(Device):
    description = 'Analog Input'
    
    @set_passed_properties(property_names = {})
    def __init__(self,name,parent_device,connection,scale_factor=1.0,units='Volts',**kwargs):
                
        self.acquisitions = []
        self.scale_factor = scale_factor
        self.units=units
        Device.__init__(self,name,parent_device,connection, **kwargs)
   
    def acquire(self,label,start_time,end_time,wait_label='',scale_factor=None,units=None):
        if scale_factor is None:
            scale_factor = self.scale_factor
        if units is None:
            units = self.units
        self.acquisitions.append({'start_time': start_time, 'end_time': end_time,
                                 'label': label, 'wait_label':wait_label, 'scale_factor':scale_factor,'units':units})
        return end_time - start_time


class Shutter(DigitalOut):
    description = 'shutter'
    
    @set_passed_properties(
        property_names = {"connection_table_properties": ["open_state"]}
        )
    def __init__(self,name,parent_device,connection,delay=(0,0),open_state=1,
                 **kwargs):

        DigitalOut.__init__(self, name, parent_device, connection, inverted=not bool(open_state), **kwargs)
        self.open_delay, self.close_delay = delay
        self.open_state = open_state
        if self.open_state == 1:
            self.allowed_states = {0: 'closed', 1: 'open'}
        elif self.open_state == 0:
            self.allowed_states = {1: 'closed', 0: 'open'}
        else:
            raise LabscriptError("Shutter %s wasn't instantiated with open_state = 0 or 1." % self.name)
        self.actual_times = {}

    # If a shutter is asked to do something at t=0, it cannot start moving
    # earlier than that.  So initial shutter states will have imprecise
    # timing. Not throwing a warning here because if I did, every run
    # would throw a warning for every shutter. The documentation will
    # have to make a point of this.
    def open(self, t):
        t_calc = t-self.open_delay if t >= self.open_delay else 0
        self.actual_times[t] = {'time': t_calc, 'instruction': 1}
        self.enable(t_calc)

    def close(self, t):
        t_calc = t-self.close_delay if t >= self.close_delay else 0
        self.actual_times[t] = {'time': t_calc, 'instruction': 0}
        self.disable(t_calc)

    def generate_code(self, hdf5_file):
        classname = self.__class__.__name__
        calibration_table_dtypes = [('name','a256'), ('open_delay',float), ('close_delay',float)]
        if classname not in hdf5_file['calibrations']:
            hdf5_file['calibrations'].create_dataset(classname, (0,), dtype=calibration_table_dtypes, maxshape=(None,))
        metadata = (self.name,self.open_delay,self.close_delay)
        dataset = hdf5_file['calibrations'][classname]
        dataset.resize((len(dataset)+1,))
        dataset[len(dataset)-1] = metadata

    def get_change_times(self, *args, **kwargs):
        retval = DigitalOut.get_change_times(self, *args, **kwargs)

        if len(self.actual_times)>1:
            sorted_times = list(self.actual_times.keys())
            sorted_times.sort()
            for i in range(len(sorted_times)-1):
                time = sorted_times[i]
                next_time = sorted_times[i+1]
                # only look at instructions that contain a state change
                if self.actual_times[time]['instruction'] != self.actual_times[next_time]['instruction']:
                    state1 = 'open' if self.actual_times[next_time]['instruction'] == 1 else 'close'
                    state2 = 'opened' if self.actual_times[time]['instruction'] == 1 else 'closed'
                    if self.actual_times[next_time]['time'] < self.actual_times[time]['time']:
                        message = "WARNING: The shutter '{:s}' is requested to {:s} too early (taking delay into account) at t={:.10f}s when it is still not {:s} from an earlier instruction at t={:.10f}s".format(self.name, state1, next_time, state2, time)
                        sys.stderr.write(message+'\n')
                elif not config.suppress_mild_warnings and not config.suppress_all_warnings:
                    state1 = 'open' if self.actual_times[next_time]['instruction'] == 1 else 'close'
                    state2 = 'opened' if self.actual_times[time]['instruction'] == 0 else 'closed'
                    message = "WARNING: The shutter '{:s}' is requested to {:s} at t={:.10f}s but was never {:s} after an earlier instruction at t={:.10f}s".format(self.name, state1, next_time, state2, time)
                    sys.stderr.write(message+'\n')
        return retval


class Trigger(DigitalOut):
    description = 'trigger device'
    allowed_states = {1:'high', 0:'low'}
    allowed_children = [TriggerableDevice]

    @set_passed_properties(property_names = {})
    def __init__(self, name, parent_device, connection, trigger_edge_type='rising',
                 **kwargs):

        DigitalOut.__init__(self,name,parent_device,connection, **kwargs)
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
        # A list of the times this trigger has been asked to trigger:
        self.triggerings = []
        
        
    def trigger(self, t, duration):
        assert duration > 0, "Negative or zero trigger duration given"
        if t != self.t0 and self.t0 not in self.instructions:
            self.disable(self.t0)
        
        start = t
        end = t + duration
        for other_start, other_duration in self.triggerings:
            other_end = other_start + other_duration
            # Check for overlapping exposures:
            if not (end < other_start or start > other_end):
                raise LabscriptError('%s %s has two overlapping triggerings: ' %(self.description, self.name) + \
                                     'one at t = %fs for %fs, and another at t = %fs for %fs.'%(start, duration, other_start, other_duration))
        self.enable(t)
        self.disable(round(t + duration,10))
        self.triggerings.append((t, duration))

    def add_device(self, device):
        if not device.connection == 'trigger':
            raise LabscriptError('The \'connection\' string of device %s '%device.name + 
                                 'to %s must be \'trigger\', not \'%s\''%(self.name, repr(device.connection)))
        DigitalOut.add_device(self, device)

        
class WaitMonitor(Trigger):
    
    @set_passed_properties(property_names={})
    def __init__(
        self,
        name,
        parent_device,
        connection,
        acquisition_device,
        acquisition_connection,
        timeout_device=None,
        timeout_connection=None,
        timeout_trigger_type='rising',
        **kwargs
    ):
        """Create a wait monitor. This is a device or devices, one of which:

        a) outputs pulses every time the master pseudoclock begins running (either at
           the start of the shot or after a wait

        b) measures the time in between those pulses in order to determine how long the
           experiment was paused for during waits

        c) optionally, produces pulses in software time that can be used to retrigger
           the master pseudoclock if a wait lasts longer than its specified timeout

        Args:

            parent_device (Device)
                The device with buffered digital outputs that should be used to produce
                the wait monitor pulses. This device must be one which is clocked by
                the master pseudoclock.

            connection (str)
                The name of the output connection of parent_device that should be used
                to produce the pulses.

            acquisition_device (Device)
                The device which is to receive those pulses as input, and that will
                measure how long between them. This does not need to be the same device
                as the wait monitor output device (corresponding to `parent_device` and
                `connection`). At time of writing, the only devices in labscript that
                can be a wait monitor acquisition device are NI DAQmx devices that have
                counter inputs.

            acquisition_connection (str)
                The name of the input connection on `acquisition_device` that is to read
                the wait monitor pulses. The user must manually connect the output
                device (`parent_device`/`connection`) to this input with a cable, in
                order that the pulses can be read by the device. On NI DAQmx devices,
                the acquisition_connection should be the name of the counter to be used
                such as 'Ctr0'. The physical connection should be made to the input
                terminal corresponding to the gate of that counter.

            timeout_device (Device, optional)
                The device that should be used to produce pulses in software time if a
                wait lasts longer than its prescribed timeout. These pulses can
                connected to the trigger input of the master pseudoclock, via a digital
                logic to 'AND' it with other trigger sources, in order to resume the
                master pseudoclock upon a wait timing out. To produce these pulses
                during a shot requires cooperation between the acquisition device and
                the timeout device code, and at present this means only NI DAQmx devices
                can be used as the timeout device (though it need not be the same device
                as the acquisition device). If not set, timeout pulses will not be
                produced and the user must manually resume the master pseudoclock via
                other means, or abort a shot if the indended resumption mechanism fails.

            timeout_connection (str, optional)
                Which connection on the timeout device should be used to produce timeout
                pulses. Since only NI DAQmx devices are supported at the moment, this
                must be a digital output on a port on the NI DAQmx device that is not
                being used. Most NI DAQmx devices have both buffered and unbuffered
                ports, so typically one would use one line of one of the unbuffered
                ports for the timeout output.

            timeout_trigger_type (str), default `'rising'`
                The edge type to be used for the timeout signal, either `'rising'` or
                `'falling'`
        """
        if compiler.wait_monitor is not None:
            raise LabscriptError("Cannot instantiate a second WaitMonitor: there can be only be one in the experiment")
        compiler.wait_monitor = self
        Trigger.__init__(self, name, parent_device, connection, **kwargs)
        if not parent_device.pseudoclock_device.is_master_pseudoclock:
            raise LabscriptError('The output device for monitoring wait durations must be clocked by the master pseudoclock device')
        
        if (timeout_device is not None) != (timeout_connection is not None):
            raise LabscriptError('Must specify both timeout_device and timeout_connection, or neither')
        self.acquisition_device = acquisition_device
        self.acquisition_connection = acquisition_connection 
        self.timeout_device = timeout_device
        self.timeout_connection = timeout_connection
        self.timeout_trigger_type = timeout_trigger_type
        
        
class DDSQuantity(Device):
    description = 'DDS'
    allowed_children = [AnalogQuantity,DigitalOut,DigitalQuantity] # Adds its own children when initialised

    @set_passed_properties(property_names = {})
    def __init__(self, name, parent_device, connection, digital_gate={}, freq_limits=None, freq_conv_class=None, freq_conv_params={},
                 amp_limits=None, amp_conv_class=None, amp_conv_params={}, phase_limits=None, phase_conv_class=None, phase_conv_params = {},
                 call_parents_add_device = True, **kwargs):
        #self.clock_type = parent_device.clock_type # Don't see that this is needed anymore
        
        # Here we set call_parents_add_device=False so that we
        # can do additional initialisation before manually calling
        # self.parent_device.add_device(self). This allows the parent's
        # add_device method to perform checks based on the code below,
        # whilst still providing us with the checks and attributes that
        # Device.__init__ gives us in the meantime.
        Device.__init__(self, name, parent_device, connection, call_parents_add_device=False, **kwargs)
                
        # Ask the parent device if it has default unit conversion classes it would like us to use:
        if hasattr(parent_device, 'get_default_unit_conversion_classes'):
            classes = self.parent_device.get_default_unit_conversion_classes(self)
            default_freq_conv, default_amp_conv, default_phase_conv = classes
            # If the user has not overridden, use these defaults. If
            # the parent does not have a default for one or more of amp,
            # freq or phase, it should return None for them.
            if freq_conv_class is None:
                freq_conv_class = default_freq_conv
            if amp_conv_class is None:
                amp_conv_class = default_amp_conv
            if phase_conv_class is None:
                phase_conv_class = default_phase_conv
        
        self.frequency = AnalogQuantity(self.name + '_freq', self, 'freq', freq_limits, freq_conv_class, freq_conv_params)
        self.amplitude = AnalogQuantity(self.name + '_amp', self, 'amp', amp_limits, amp_conv_class, amp_conv_params)
        self.phase = AnalogQuantity(self.name + '_phase', self, 'phase', phase_limits, phase_conv_class, phase_conv_params)

        self.gate = None
        if 'device' in digital_gate and 'connection' in digital_gate:            
            self.gate = DigitalOut(name + '_gate', digital_gate['device'], digital_gate['connection'])
        # Did they only put one key in the dictionary, or use the wrong keywords?
        elif len(digital_gate) > 0:
            raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.' % (self.name))
        
        # If the user has not specified a gate, and the parent device
        # supports gating of DDS output, it should add a gate to this
        # instance in its add_device method, which is called below. If
        # they *have* specified a gate device, but the parent device
        # has its own gating (such as the PulseBlaster), it should
        # check this and throw an error in its add_device method. See
        # labscript_devices.PulseBlaster.PulseBlaster.add_device for an
        # example of this.
        # In some subclasses we need to hold off on calling the parent
        # device's add_device function until further code has run,
        # e.g., see PulseBlasterDDS in PulseBlaster.py
        if call_parents_add_device:
            self.parent_device.add_device(self)
        
    def setamp(self, t, value, units=None):
        self.amplitude.constant(t, value, units)
        
    def setfreq(self, t, value, units=None):
        self.frequency.constant(t, value, units)
        
    def setphase(self, t, value, units=None):
        self.phase.constant(t, value, units)
        
    def enable(self, t):
        if self.gate is None:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.' % (self.name))
        self.gate.go_high(t)

    def disable(self, t):
        if self.gate is None:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.' % (self.name))
        self.gate.go_low(t)
            
    def pulse(self, t, duration, amplitude, frequency, phase=None, amplitude_units = None, frequency_units = None, phase_units = None, print_summary=False):
        if print_summary:
            functions.print_time(t, '%s pulse at %.4f MHz for %.3f ms' % (self.name, frequency/MHz, duration/ms))
        self.setamp(t, amplitude, amplitude_units)
        if frequency is not None:
            self.setfreq(t, frequency, frequency_units)
        if phase is not None:
            self.setphase(t, phase, phase_units)
        if amplitude != 0 and self.gate is not None:
            self.enable(t)
            self.disable(t + duration)
            self.setamp(t + duration, 0)
        return duration

class DDS(DDSQuantity):
    pass

class StaticDDS(Device):
    description = 'Static RF'
    allowed_children = [StaticAnalogQuantity,DigitalOut,StaticDigitalOut]
    
    @set_passed_properties(property_names = {})
    def __init__(self,name,parent_device,connection,digital_gate = {},freq_limits = None,freq_conv_class = None,freq_conv_params = {},amp_limits=None,amp_conv_class = None,amp_conv_params = {},phase_limits=None,phase_conv_class = None,phase_conv_params = {},
                 **kwargs):
        #self.clock_type = parent_device.clock_type # Don't see that this is needed anymore
        
        # We tell Device.__init__ to not call
        # self.parent.add_device(self), we'll do that ourselves later
        # after further intitialisation, so that the parent can see the
        # freq/amp/phase objects and manipulate or check them from within
        # its add_device method.
        Device.__init__(self,name,parent_device,connection, call_parents_add_device=False, **kwargs)

        # Ask the parent device if it has default unit conversion classes it would like us to use:
        if hasattr(parent_device, 'get_default_unit_conversion_classes'):
            classes = parent_device.get_default_unit_conversion_classes(self)
            default_freq_conv, default_amp_conv, default_phase_conv = classes
            # If the user has not overridden, use these defaults. If
            # the parent does not have a default for one or more of amp,
            # freq or phase, it should return None for them.
            if freq_conv_class is None:
                freq_conv_class = default_freq_conv
            if amp_conv_class is None:
                amp_conv_class = default_amp_conv
            if phase_conv_class is None:
                phase_conv_class = default_phase_conv

        self.frequency = StaticAnalogQuantity(self.name+'_freq',self,'freq',freq_limits,freq_conv_class,freq_conv_params)
        self.amplitude = StaticAnalogQuantity(self.name+'_amp',self,'amp',amp_limits,amp_conv_class,amp_conv_params)
        self.phase = StaticAnalogQuantity(self.name+'_phase',self,'phase',phase_limits,phase_conv_class,phase_conv_params)        
        
        if 'device' in digital_gate and 'connection' in digital_gate:            
            self.gate = DigitalOut(self.name+'_gate',digital_gate['device'],digital_gate['connection'])
        # Did they only put one key in the dictionary, or use the wrong keywords?
        elif len(digital_gate) > 0:
            raise LabscriptError('You must specify the "device" and "connection" for the digital gate of %s.'%(self.name))
        # Now we call the parent's add_device method. This is a must, since we didn't do so earlier from Device.__init__.
        self.parent_device.add_device(self)
        
    def setamp(self,value,units=None):
        self.amplitude.constant(value,units)
        
    def setfreq(self,value,units=None):
        self.frequency.constant(value,units)
        
    def setphase(self,value,units=None):
        self.phase.constant(value,units) 
            
    def enable(self,t=None):        
        if self.gate:
            self.gate.go_high(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the enable(t) method.'%(self.name))
                        
    def disable(self,t=None):
        if self.gate:
            self.gate.go_low(t)
        else:
            raise LabscriptError('DDS %s does not have a digital gate, so you cannot use the disable(t) method.'%(self.name))
              
class LabscriptError(Exception):
    pass

def save_time_markers(hdf5_file):
    time_markers = compiler.time_markers
    dtypes = [('label','a256'), ('time', float), ('color', '(1,3)int')]
    data_array = zeros(len(time_markers), dtype=dtypes)
    for i, t in enumerate(time_markers):
        data_array[i] = time_markers[t]["label"], t, time_markers[t]["color"]
    time_markers_dataset = hdf5_file.create_dataset('time_markers', data = data_array)

def generate_connection_table(hdf5_file):
    connection_table = []
    devicedict = {}
    
    # Only use a string dtype as long as is needed:
    max_BLACS_conn_length = -1

    for device in compiler.inventory:
        devicedict[device.name] = device

        unit_conversion_parameters = device._properties['unit_conversion_parameters']
        serialised_unit_conversion_parameters = labscript_utils.properties.serialise(unit_conversion_parameters)

        properties = device._properties["connection_table_properties"]
        serialised_properties = labscript_utils.properties.serialise(properties)
        
        # If the device has a BLACS_connection atribute, then check to see if it is longer than the size of the hdf5 column
        if hasattr(device,"BLACS_connection"):
            # Make sure it is a string!
            BLACS_connection = str(device.BLACS_connection)
            if len(BLACS_connection) > max_BLACS_conn_length:
                max_BLACS_conn_length = len(BLACS_connection)
        else:
            BLACS_connection = ""
            
            #if there is no BLACS connection, make sure there is no "gui" or "worker" entry in the connection table properties
            if 'worker' in properties or 'gui' in properties:
                raise LabscriptError('You cannot specify a remote GUI or worker for a device (%s) that does not have a tab in BLACS'%(device.name))
            
            
        connection_table.append((device.name, device.__class__.__name__,
                                 device.parent_device.name if device.parent_device else str(None),
                                 str(device.connection if device.parent_device else str(None)),
                                 device.unit_conversion_class.__name__ if hasattr(device,"unit_conversion_class") and device.unit_conversion_class is not None else str(None),
                                 serialised_unit_conversion_parameters,
                                 BLACS_connection,
                                 serialised_properties))
    
    connection_table.sort()
    vlenstring = h5py.special_dtype(vlen=str)
    connection_table_dtypes = [('name','a256'), ('class','a256'), ('parent','a256'), ('parent port','a256'),
                               ('unit conversion class','a256'), ('unit conversion params', vlenstring),
                               ('BLACS_connection','a'+str(max_BLACS_conn_length)),
                               ('properties', vlenstring)]
    connection_table_array = empty(len(connection_table),dtype=connection_table_dtypes)
    for i, row in enumerate(connection_table):
        connection_table_array[i] = row
    dataset = hdf5_file.create_dataset('connection table', compression=config.compression, data=connection_table_array, maxshape=(None,))
    
    if compiler.master_pseudoclock is None:
        master_pseudoclock_name = 'None'
    else:
        master_pseudoclock_name = compiler.master_pseudoclock.name
    dataset.attrs['master_pseudoclock'] = master_pseudoclock_name
  
  
def save_labscripts(hdf5_file):
    if compiler.labscript_file is not None:
        script_text = open(compiler.labscript_file).read()
    else:
        script_text = ''
    script = hdf5_file.create_dataset('script',data=script_text)
    script.attrs['name'] = os.path.basename(compiler.labscript_file).encode() if compiler.labscript_file is not None else ''
    script.attrs['path'] = os.path.dirname(compiler.labscript_file).encode() if compiler.labscript_file is not None else sys.path[0]
    try:
        import labscriptlib
        prefix = os.path.dirname(labscriptlib.__file__)
        for module in sys.modules.values():
            if getattr(module, '__file__', None) is not None:
                path = os.path.abspath(module.__file__)
                if path.startswith(prefix) and (path.endswith('.pyc') or path.endswith('.py')):
                    path = path.replace('.pyc', '.py')
                    save_path = 'labscriptlib/' + path.replace(prefix, '').replace('\\', '/').replace('//', '/')
                    if save_path in hdf5_file:
                        # Don't try to save the same module script twice! 
                        # (seems to at least double count __init__.py when you import an entire module as in from labscriptlib.stages import * where stages is a folder with an __init__.py file.
                        # Doesn't seem to want to double count files if you just import the contents of a file within a module
                        continue
                    hdf5_file.create_dataset(save_path, data=open(path).read())
                    if compiler.save_hg_info:
                        hg_commands = [['log', '--limit', '1'], ['status'], ['diff']]
                        for command in hg_commands:
                            process = subprocess.Popen(['hg'] + command + [os.path.split(path)[1]], cwd=os.path.split(path)[0],
                                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, startupinfo=startupinfo)
                            info, err = process.communicate()
                            if info or err:
                                hdf5_file[save_path].attrs['hg ' + str(command[0])] = info.decode('utf-8') + '\n' + err.decode('utf-8')
    except ImportError:
        pass
    except WindowsError if os.name == 'nt' else None:
        sys.stderr.write('Warning: Cannot save Mercurial data for imported scripts. Check that the hg command can be run from the command line.\n')


def write_device_properties(hdf5_file):
    for device in compiler.inventory:
        device_properties = device._properties["device_properties"]

        # turn start_order and stop_order into device_properties,
        # but only if the device has a BLACS_connection attribute. start_order
        # or stop_order being None means the default order, zero. An order
        # being specified on a device without a BLACS_connection is an error.
        for attr_name in ['start_order', 'stop_order']:
            attr = getattr(device, attr_name)
            if hasattr(device, 'BLACS_connection'):
                if attr is None:
                    # Default:
                    attr = 0
                device_properties[attr_name] = attr
            elif attr is not None:
                msg = ('cannot set %s on device %s, which does ' % (attr_name, device.name) +
                       'not have a BLACS_connection attribute and thus is not started and stopped directly by BLACS.')
                raise TypeError(msg)

        # Special case: We don't create the group if the only property is an
        # empty dict called 'added properties'. This is because this property
        # is present in all devices, and represents a place to pass in
        # arbitrary data from labscript experiment scripts. We don't want a
        # group for every device if nothing is actually being passed in, so we
        # ignore this case.
        if device_properties and device_properties != {'added_properties':{}}:
            # Create group if doesn't exist:
            if not device.name in hdf5_file['devices']:
                hdf5_file['/devices'].create_group(device.name)
            labscript_utils.properties.set_device_properties(hdf5_file, device.name, device_properties)


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
        if compiler.wait_monitor.timeout_device is None:
            timeout_device = ''
        else:
            timeout_device = compiler.wait_monitor.timeout_device.name
        if compiler.wait_monitor.timeout_connection is None:
            timeout_connection = ''
        else:
            timeout_connection = compiler.wait_monitor.timeout_connection
        timeout_trigger_type = compiler.wait_monitor.timeout_trigger_type
    else:
        acquisition_device, acquisition_connection, timeout_device, timeout_connection, timeout_trigger_type = '', '', '', '', ''
    dataset.attrs['wait_monitor_acquisition_device'] = acquisition_device
    dataset.attrs['wait_monitor_acquisition_connection'] = acquisition_connection
    dataset.attrs['wait_monitor_timeout_device'] = timeout_device
    dataset.attrs['wait_monitor_timeout_connection'] = timeout_connection
    dataset.attrs['wait_monitor_timeout_trigger_type'] = timeout_trigger_type

    
def generate_code():
    if compiler.hdf5_filename is None:
        raise LabscriptError('hdf5 file for compilation not set. Please call labscript_init')
    elif not os.path.exists(compiler.hdf5_filename):
        with h5py.File(compiler.hdf5_filename ,'w') as hdf5_file:
            hdf5_file.create_group('globals')
    with h5py.File(compiler.hdf5_filename, 'a') as hdf5_file:
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
            if device.parent_device is None:
                device.generate_code(hdf5_file)
                
        save_time_markers(hdf5_file)
        generate_connection_table(hdf5_file)
        write_device_properties(hdf5_file)
        generate_wait_table(hdf5_file)
        save_labscripts(hdf5_file)

        # Save shot properties:
        group = hdf5_file.create_group('shot_properties')
        set_attributes(group, compiler.shot_properties)


def trigger_all_pseudoclocks(t='initial'):
    # Must wait this long before providing a trigger, in case child clocks aren't ready yet:
    wait_delay = compiler.wait_delay
    if type(t) in [str, bytes] and t == 'initial':
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
    # So find the max of 1.0/clock_limit of every clockline on every pseudoclock of the master pseudoclock
    master_pseudoclock_delay = max(1.0/compiler.master_pseudoclock.clock_limit, max_or_zero([1.0/clockline.clock_limit for pseudoclock in compiler.master_pseudoclock.child_devices for clockline in pseudoclock.child_devices]))
    max_delay = max(compiler.trigger_duration + master_pseudoclock_delay, max_delay_time)    
    return max_delay + wait_delay
    
def wait(label, t, timeout=5):
    if not str(label):
        raise LabscriptError('Wait must have a name')
    max_delay = trigger_all_pseudoclocks(t)
    if t in compiler.wait_table:
        raise LabscriptError('There is already a wait at t=%s'%str(t))
    if any([label==existing_label for existing_label, _ in compiler.wait_table.values()]):
        raise LabscriptError('There is already a wait named %s'%str(label))
    compiler.wait_table[t] = str(label), float(timeout)
    return max_delay

def add_time_marker(t, label, color=None, verbose=False):
    """Add a marker for the specified time. These markers are saved in the HDF5 file.
    This allows one to label that time with a string label, and a color that
    applications may use to represent this part of the experiment. The color may be
    specified as an RGB tuple, or a string of the color name such as 'red', or a string
    of its hex value such as '#ff88g0'. If verbose=True, this funtion also prints the
    label and time, which can be useful to view during shot compilation.

    Runviewer displays these markers and allows one to manipulate the time axis based on
    them, and BLACS' progress bar plugin displays the name and colour of the most recent
    marker as the shot is running"""
    if isinstance(color, (str, bytes)):
        import PIL.ImageColor
        color = PIL.ImageColor.getrgb(color)
    if color is None:
        color = (-1, -1, -1)
    elif not all(0 <= n <= 255 for n in color):
        raise ValueError("Invalid RGB tuple %s" % str(color))
    if verbose:
        functions.print_time(t, label)
    compiler.time_markers[t] = {"label": label, "color": color}

def start():
    compiler.start_called = True
    # Get and save some timing info about the pseudoclocks:
    # TODO: What if you need to trigger individual Pseudolocks on the one device, rather than the PseudoclockDevice as a whole?
    pseudoclocks = [device for device in compiler.inventory if isinstance(device, PseudoclockDevice)]
    compiler.all_pseudoclocks = pseudoclocks
    toplevel_devices = [device for device in compiler.inventory if device.parent_device is None]
    master_pseudoclocks = [pseudoclock for pseudoclock in pseudoclocks if pseudoclock.is_master_pseudoclock]
    if len(master_pseudoclocks) > 1:
        raise LabscriptError('Cannot have more than one master pseudoclock')
    if not toplevel_devices:
        raise LabscriptError('No toplevel devices and no master pseudoclock found')
    elif pseudoclocks:
        (master_pseudoclock,) = master_pseudoclocks
        compiler.master_pseudoclock = master_pseudoclock
        # Which pseudoclock requires the longest pulse in order to trigger it?
        compiler.trigger_duration = max_or_zero([pseudoclock.trigger_minimum_duration for pseudoclock in pseudoclocks if not pseudoclock.is_master_pseudoclock])
        
        trigger_clock_limits = [pseudoclock.trigger_device.clock_limit for pseudoclock in pseudoclocks if not pseudoclock.is_master_pseudoclock]
        if len(trigger_clock_limits) > 0:
            min_clock_limit = min(trigger_clock_limits)
            min_clock_limit = min([min_clock_limit, master_pseudoclock.clock_limit])
        else:
            min_clock_limit = master_pseudoclock.clock_limit
    
        # ensure we don't tick faster than attached devices to the master pseudoclock can handle
        clockline_limits = [clockline.clock_limit for pseudoclock in master_pseudoclock.child_devices for clockline in pseudoclock.child_devices]
        min_clock_limit = min(min_clock_limit, min(clockline_limits))
    
        # check the minimum trigger duration for the waitmonitor
        if compiler.wait_monitor is not None:
            wait_monitor_minimum_pulse_width = getattr(
                compiler.wait_monitor.acquisition_device,
                'wait_monitor_minimum_pulse_width',
                0,
            )
            compiler.trigger_duration = max(
                compiler.trigger_duration, wait_monitor_minimum_pulse_width
            )
            
        # Provide this, or the minimum possible pulse, whichever is longer:
        compiler.trigger_duration = max(2.0/min_clock_limit, compiler.trigger_duration) + 2*master_pseudoclock.clock_resolution
        # Must wait this long before providing a trigger, in case child clocks aren't ready yet:
        compiler.wait_delay = max_or_zero([pseudoclock.wait_delay for pseudoclock in pseudoclocks if not pseudoclock.is_master_pseudoclock])
        
        # Have the master clock trigger pseudoclocks at t = 0:
        max_delay = trigger_all_pseudoclocks()
    else:
        # No pseudoclocks, only other toplevel devices:
        compiler.master_pseudoclock = None
        compiler.trigger_duration = 0
        compiler.wait_delay = 0
        max_delay = 0
    return max_delay
    
def stop(t, target_cycle_time=None, cycle_time_delay_after_programming=False):
    """Indicate the end of an experiment at the given time, and initiate compilation of
    instructions, saving them to the HDF5 file. Configure some shot options.
    
    t (float or None), default: `None`
        The end time of the experiment.

    target_cycle_time (float or None), default: `None`
        How long, in seconds, after the previous shot was started, should this shot be
        started by BLACS. This allows one to run shots at a constant rate even if they
        are of different durations. If `None`, BLACS will run the next shot immediately
        after the previous shot completes. Otherwise, BLACS will delay starting this
        shot until the cycle time has elapsed. This is a request only, and may not be
        met if running/programming/saving data from a shot takes long enough that it
        cannot be met. This functionality requires the BLACS `cycle_time` plugin to be
        enabled in labconfig. Its accuracy is also limited by software timing,
        requirements of exact cycle times beyond software timing should be instead done
        using hardware triggers to Pseudoclocks.

    cycle_time_delay_after_programming (bool), default: `False`
        Whether the BLACS cycle_time plugin should insert the required delay for the
        target cycle time *after* programming devices, as opposed to before programming
        them. This is more precise, but may cause some devices to output their first
        instruction for longer than desired, since some devices begin outputting their
        first instruction as soon as they are programmed rather than when they receive
        their first clock tick. If not set, the *average* cycle time will still be just
        as close to as requested (so long as there is adequate time available), however
        the time interval between the same part of the experiment from one shot to the
        next will not be as precise due to variations in programming time.
    """
    # Indicate the end of an experiment and initiate compilation:
    if t == 0:
        raise LabscriptError('Stop time cannot be t=0. Please make your run a finite duration')
    for device in compiler.inventory:
        if isinstance(device, PseudoclockDevice):
            device.stop_time = t
    if target_cycle_time is not None:
        # Ensure we have a valid type for this
        target_cycle_time = float(target_cycle_time)
    compiler.shot_properties['target_cycle_time'] = target_cycle_time
    compiler.shot_properties['cycle_time_delay_after_programming'] = cycle_time_delay_after_programming
    generate_code()

# TO_DELETE:runmanager-batchompiler-agnostic
#   entire function load_globals can be deleted
def load_globals(hdf5_filename):
    import runmanager
    params = runmanager.get_shot_globals(hdf5_filename)
    with h5py.File(hdf5_filename,'r') as hdf5_file:
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
            # 'None' is stored as an h5py null object reference:
            if isinstance(params[name], h5py.Reference) and not params[name]:
                params[name] = None
            _builtins_dict[name] = params[name]
            
# TO_DELETE:runmanager-batchompiler-agnostic 
#   load_globals_values=True            
def labscript_init(hdf5_filename, labscript_file=None, new=False, overwrite=False, load_globals_values=True):
    # save the builtins for later restoration in labscript_cleanup
    compiler._existing_builtins_dict = _builtins_dict.copy()
    
    if new:
        # defer file creation until generate_code(), so that filesystem
        # is not littered with h5 files when the user merely imports
        # labscript. If the file already exists, and overwrite is true, delete it so we get one fresh.
        if os.path.exists(hdf5_filename) and overwrite:
            os.unlink(hdf5_filename)
    elif not os.path.exists(hdf5_filename):
        raise LabscriptError('Provided hdf5 filename %s doesn\'t exist.'%hdf5_filename)
    # TO_DELETE:runmanager-batchompiler-agnostic 
    elif load_globals_values:
        load_globals(hdf5_filename) 
    # END_DELETE:runmanager-batchompiler-agnostic 
    
    compiler.hdf5_filename = hdf5_filename
    if labscript_file is None:
        import __main__
        labscript_file = __main__.__file__
    compiler.labscript_file = os.path.abspath(labscript_file)
    

def labscript_cleanup():
    """restores builtins and the labscript module to its state before
    labscript_init() was called"""
    for name in _builtins_dict.copy(): 
        if name not in compiler._existing_builtins_dict:
            del _builtins_dict[name]
        else:
            _builtins_dict[name] = compiler._existing_builtins_dict[name]
    
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
    compiler.time_markers = {}
    compiler._PrimaryBLACS = None
    compiler.save_hg_info = True
    compiler.shot_properties = {}

class compiler(object):
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
    time_markers = {}
    _PrimaryBLACS = None
    save_hg_info = True
    shot_properties = {}

    # safety measure in case cleanup is called before init
    _existing_builtins_dict = _builtins_dict.copy() 
