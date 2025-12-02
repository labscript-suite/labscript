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

"""Everything else including the `start()`, `stop()`, and `wait()` functions - all other
classes are also imported here for backwards compatibility"""

import builtins
import os
import sys
import subprocess
import keyword
import threading
from functools import lru_cache
import numpy as np

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
from labscript_utils.filewatcher import FileWatcher

# This imports the default Qt library that other labscript suite code will
# import as well, since it all uses qtutils. By having a Qt library already
# imported, we ensure matplotlib (imported by pylab) will notice this and use
# the same Qt library and API version, and therefore not conflict with any
# other code is using:
import qtutils

from pylab import *

from labscript_utils import dedent
from labscript_utils.properties import set_attributes

import labscript.functions as functions
try:
    from labscript_utils.unitconversions import *
except ImportError:
    sys.stderr.write('Warning: Failed to import unit conversion classes\n')

# Imports to maintain backwards compatibility with previous module interface
from .base import Device
from .compiler import compiler
from .constants import *
from .core import (
    ClockLine,
    IntermediateDevice,
    Pseudoclock,
    PseudoclockDevice,
    TriggerableDevice,
)
from .inputs import AnalogIn
from .outputs import (
    Output,
    AnalogOut,
    AnalogQuantity,
    DDS,
    DDSQuantity,
    DigitalOut,
    DigitalQuantity,
    Shutter,
    StaticDDS,
    StaticAnalogOut,
    StaticAnalogQuantity,
    StaticDigitalOut,
    StaticDigitalQuantity,
    Trigger,
)
from .utils import (
    LabscriptError,
    bitfield,
    fastflatten,
    max_or_zero,
    set_passed_properties,
    suppress_all_warnings,
    suppress_mild_warnings
)

# Create a reference to the builtins dict 
# update this if accessing builtins ever changes
_builtins_dict = builtins.__dict__
    
# Startupinfo, for ensuring subprocesses don't launch with a visible command window:
if os.name=='nt':
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= 1 #subprocess.STARTF_USESHOWWINDOW # This variable isn't defined, but apparently it's equal to one.
else:
    startupinfo = None

# Alias for backwards compatibility. The config has been moved into the compiler with
# the other settings.
config = compiler 
no_warnings = suppress_all_warnings # Historical alias


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
        """Create a wait monitor. 

        This is a device or devices, one of which:

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


def save_time_markers(hdf5_file):
    """Save shot time markers to the shot file.

    Args:
        hdf5_file (:obj:`h5py:h5py.File`): Handle to file to save to.
    """
    time_markers = compiler.time_markers
    dtypes = [('label','a256'), ('time', float), ('color', '(1,3)int')]
    data_array = zeros(len(time_markers), dtype=dtypes)
    for i, t in enumerate(time_markers):
        data_array[i] = time_markers[t]["label"], t, time_markers[t]["color"]
    time_markers_dataset = hdf5_file.create_dataset('time_markers', data = data_array)

def generate_connection_table(hdf5_file):
    """Generates the connection table for the compiled shot.

    Args:
        hdf5_file (:obj:`h5py:h5py.File`): Handle to file to save to.
    """
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
          
        if getattr(device, 'unit_conversion_class', None) is not None:
            c = device.unit_conversion_class
            unit_conversion_class_repr = f"{c.__module__}.{c.__name__}"
        else:
            unit_conversion_class_repr = repr(None)

        connection_table.append((device.name, device.__class__.__name__,
                                 device.parent_device.name if device.parent_device else str(None),
                                 str(device.connection if device.parent_device else str(None)),
                                 unit_conversion_class_repr,
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
    dataset = hdf5_file.create_dataset('connection table', compression=compiler.compression, data=connection_table_array, maxshape=(None,))
    
    if compiler.master_pseudoclock is None:
        master_pseudoclock_name = 'None'
    else:
        master_pseudoclock_name = compiler.master_pseudoclock.name
    dataset.attrs['master_pseudoclock'] = master_pseudoclock_name

# Create a dictionary for caching results from vcs commands. The keys will be
# the paths to files that are saved during save_labscripts(). The values will be
# a list of tuples of the form (command, info, err); see the "Returns" section
# of the _run_vcs_commands() docstring for more info. Also create a FileWatcher
# instance for tracking when vcs results need updating. The callback will
# replace the outdated cache entry with a new list of updated vcs commands and
# outputs.
_vcs_cache = {}
_vcs_cache_rlock = threading.RLock()
def _file_watcher_callback(name, info, event):
    with _vcs_cache_rlock:
        _vcs_cache[name] = _run_vcs_commands(name)

_file_watcher = FileWatcher(_file_watcher_callback)

@lru_cache(None)
def _have_vcs(vcs):
    try:
        subprocess.check_output([vcs, '--version'])
        return True
    except FileNotFoundError:
        msg = f"""Warning: Cannot run {vcs} commands, {vcs} info for labscriptlib files
        will not be saved. You can disable this warning by setting
        [labscript]/save_{vcs}_info = False in labconfig."""
        sys.stderr.write(dedent(msg) + '\n')
        return False

def _run_vcs_commands(path):
    """Run some VCS commands on a file and return their output.
    
    The function is used to gather up version control system information so that
    it can be stored in the hdf5 files of shots. This is for convenience and
    compliments the full copy of the file already included in the shot file.
    
    Whether hg and git commands are run is controlled by the `save_hg_info`
    and `save_git_info` options in the `[labscript]` section of the labconfig.

    Args:
        path (str): The path with file name and extension of the file on which
            the commands will be run. The working directory will be set to the
            directory containing the specified file.

    Returns:
        results (list of (tuple, str, str)): A list of tuples, each
            containing information related to one vcs command of the form
            (command, info, err). The first entry in that tuple is itself a
            tuple of strings which was passed to subprocess.Popen() in order to
            run the command. Then info is a string that contains the text
            printed to stdout by that command, and err contains the text printed
            to stderr by the command.
    """
    # Gather together a list of commands to run.
    module_directory, module_filename = os.path.split(path)
    vcs_commands = []
    if compiler.save_hg_info and _have_vcs('hg'):
        hg_commands = [
            ['log', '--limit', '1'],
            ['status'],
            ['diff'],
        ]
        for command in hg_commands:
            command = tuple(['hg'] + command + [module_filename])
            vcs_commands.append((command, module_directory))
    if compiler.save_git_info and _have_vcs('git'):
        git_commands = [
            ['branch', '--show-current'],
            ['describe', '--tags', '--always', 'HEAD'],
            ['rev-parse', 'HEAD'],
            ['diff', 'HEAD', module_filename],
        ]
        for command in git_commands:
            command = tuple(['git'] + command)
            vcs_commands.append((command, module_directory))

    # Now go through and start running the commands.
    process_list = []
    for command, module_directory in vcs_commands:
        process = subprocess.Popen(
            command,
            cwd=module_directory,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            startupinfo=startupinfo,
        )
        process_list.append((command, process))

    # Gather up results from the commands issued.
    results = []
    for command, process in process_list:
        info, err = process.communicate()
        info = info.decode('utf-8')
        err = err.decode('utf-8')
        results.append((command, info, err))
    return results
 
def save_labscripts(hdf5_file):
    """Writes the script files for the compiled shot to the shot file.

    If `save_hg_info` labconfig parameter is `True`, will attempt to save
    hg version info as an attribute.

    Args:
        hdf5_file (:obj:`h5py:h5py.File`): Handle to file to save to.
    """
    if compiler.labscript_file is not None:
        script_text = open(compiler.labscript_file).read()
    else:
        script_text = ''
    script = hdf5_file.create_dataset('script',data=script_text)
    script.attrs['name'] = os.path.basename(compiler.labscript_file).encode() if compiler.labscript_file is not None else ''
    script.attrs['path'] = os.path.dirname(compiler.labscript_file).encode() if compiler.labscript_file is not None else sys.path[0]
    try:
        import labscriptlib
    except ImportError:
        return
    prefix = os.path.dirname(labscriptlib.__file__)
    for module in sys.modules.values():
        if getattr(module, '__file__', None) is not None:
            path = os.path.abspath(module.__file__)
            if path.startswith(prefix) and (path.endswith('.pyc') or path.endswith('.py')):
                if 'signature_bootstrap.py' in path or 'shibokensupport' in path:
                    # ignore PySide6 shenaniganry
                    continue
                path = path.replace('.pyc', '.py')
                save_path = 'labscriptlib/' + path.replace(prefix, '').replace('\\', '/').replace('//', '/')
                if save_path in hdf5_file:
                    # Don't try to save the same module script twice! (seems to at least
                    # double count __init__.py when you import an entire module as in
                    # from labscriptlib.stages import * where stages is a folder with an
                    # __init__.py file. Doesn't seem to want to double count files if
                    # you just import the contents of a file within a module
                    continue
                hdf5_file.create_dataset(save_path, data=open(path).read())
                with _vcs_cache_rlock:
                    if not path in _vcs_cache:
                        # Add file to watch list and create its entry in the cache.
                        _file_watcher.add_file(path)
                        _file_watcher_callback(path, None, None)
                    # Save the cached vcs output to the file.
                    for command, info, err in _vcs_cache[path]:
                        attribute_str = command[0] + ' ' + command[1]
                        hdf5_file[save_path].attrs[attribute_str] = (info + '\n' + err)



def write_device_properties(hdf5_file):
    """Writes device_properties for each device in compiled shot to shto file.

    Args:
        hdf5_file (:obj:`h5py:h5py.File`): Handle to file to save to.
    """
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
    """Generates the wait table for the shot and saves it to the shot file.

    Args:
        hdf5_file (:obj:`h5py:h5py.File`): Handle to file to save to.
    """
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
    """Compiles a shot and saves it to the shot file.
    """
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
    """Commands pseudoclocks to pause until resumed by an external trigger, or a timeout is reached.

    Args:
        label (str): Unique name for wait.
        t (float): Time, in seconds, at which experiment should begin the wait.
        timeout (float, optional): Maximum length of the wait, in seconds. After
            this time, the pseudoclocks are automatically re-triggered.

    Returns:
        float: Time required for all pseudoclocks to resume execution once
        wait has completed.            
    """
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
    """Indicates the end of the connection table and the start of the 
    experiment logic.

    Returns:
        float: Time required for all pseudoclocks to start execution.
    """
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
    instructions, saving them to the HDF5 file. Configures some shot options.
    
    Args:
        t (float): The end time of the experiment.

        target_cycle_time (float, optional): default: `None`
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

        cycle_time_delay_after_programming (bool, optional): default: `False`
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
    import labscript_utils.shot_utils
    params = labscript_utils.shot_utils.get_shot_globals(hdf5_filename)
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
            if type(params[name]) == np.bool_:
                params[name] = bool(params[name])
            # 'None' is stored as an h5py null object reference:
            if isinstance(params[name], h5py.Reference) and not params[name]:
                params[name] = None
            _builtins_dict[name] = params[name]
            
# TO_DELETE:runmanager-batchompiler-agnostic 
#   load_globals_values=True            
def labscript_init(hdf5_filename, labscript_file=None, new=False, overwrite=False, load_globals_values=True):
    """Initialises labscript and prepares for compilation.

    Args:
        hdf5_filename (str): Path to shot file to compile.
        labscript_file: Handle to the labscript file.
        new (bool, optional): If `True`, ensure a new shot file is created.
        overwrite (bool, optional): If `True`, overwrite existing shot file, if it exists.
        load_globals_values (bool, optional): If `True`, load global values 
            from the existing shot file.
    """
    # save the builtins for later restoration in labscript_cleanup
    compiler.save_builtins_state()
    
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
    compiler.reset()
