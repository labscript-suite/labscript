#####################################################################
#                                                                   #
# /devices/adwin.py                                                 #
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
from labscript import Device, PseudoClock, IntermediateDevice, AnalogOut, DigitalOut, bitfield
from pylab import *
default_cycle_time = 2500/300e6 # 8.333us

# Notes:
# The ADWin runs at 300MHz. The cycle time should be specified in hardware programming in units of this clock speed.
# Subsequently, instruction timing must be specified in units of cycles.

# Voltages are specified with a 16 bit unsigned integer, mapping the range [-10,10) volts.
# There are 32 digital outputs on a card
# There are 8 analog outs on a card

from IPython import embed

class Instruction(dict):
    """A dictionary-like class that provides the interface to instructions
    that the rest of labscript expects, (ie, a dictionary). This is a
    temporary measure until labscript is updated to use objects for all
    instructions, rather than dictionarys."""
    
    def __getitem__(self, key):
        if key == 'initial time':
            return self.initial_time
        elif key == 'description':
            return self.__class__.__name__
        elif key == 'duration':
            return self.duration
        elif key == 'function':
            return self.evaluate
        elif key == 'clock rate':
            return self.clock_rate
        elif key == 'end time':
            return self.end_time
        elif key == 'units':
            return self.units
        else:
            import IPython
            IPython.embed()
            raise KeyError(key)

    def __setitem__(self, key, value):
        if key == 'initial time':
            self.initial_time = value
        elif key == 'duration':
            self.duration = value
        elif key == 'clock rate':
            self.clock_rate = value
        elif key == 'end time':
            self.end_time = value
        else:
            raise KeyError(key)

class RampInstruction(Instruction):
    def __init__(self, t, duration, A, B, C, clock_rate, units=None):
        self.initial_time = t
        self.duration = duration
        self.end_time = self.initial_time + duration
        self.A = A
        self.B = B
        self.C = C
        self.clock_rate = clock_rate
        self.units = units
    
    def copy(self):
        # Override dict.copy() so we get the same class back when
        # labscript copies an instruction. We can't call __init__, as
        # we don't know the call signature of subclasses. We have to
        # call __new__ and then manually set attributes.
        new_instruction = self.__class__.__new__(self.__class__)
        new_instruction.initial_time = self.initial_time
        new_instruction.duration = self.duration
        new_instruction.end_time = self.end_time
        new_instruction.A = self.A
        new_instruction.B = self.B
        new_instruction.C = self.C
        new_instruction.clock_rate = self.clock_rate
        new_instruction.units = self.units
        return new_instruction
        
    def __repr__(self):
        return object.__repr__(self)
        
    def evaluate(self,t):
        # To be overridden by subclasses
        raise NotImplementedError
        
class LinearRamp(RampInstruction):
    def __init__(self, t, duration, initial_value, final_value, clock_rate):
        RampInstruction.__init__(self, t, duration, A=final_value, B=0, C=initial_value, clock_rate=clock_rate)

    def evaluate(self, t):
        return (self.A - self.C)/self.duration * t  + self.C
        
class SinRamp(RampInstruction):
    def __init__(self, t, duration, amplitude, offset, angular_period, clock_rate):
        RampInstruction.__init__(self, t, duration, A=amplitude, B=angular_period, C=offset, clock_rate=clock_rate)

    def evaluate(self, t):
        return self.A*sin(t/self.B) + self.C
        
class CosRamp(RampInstruction):
    def __init__(self, t, duration, amplitude, offset, angular_period, clock_rate):
        RampInstruction.__init__(self, t, duration, A=amplitude, B=angular_period, C=offset, clock_rate=clock_rate)

    def evaluate(self, t):
        return self.A*cos(t/self.B) + self.C
        
class ExpRamp(RampInstruction):
    def __init__(self, t, duration, amplitude, offset, time_constant, clock_rate):
        RampInstruction.__init__(self, t, duration, A=amplitude, B=time_constant, C=offset, clock_rate=clock_rate)

    def evaluate(self, t):
        return self.A*exp(t/self.B) + self.C
        
        
class ADWinAnalogOut(AnalogOut):
    def linear_ramp(self, t, duration, initial, final):
        instruction = LinearRamp(t, duration, initial, final, self.parent_device.clock_limit)
        self.add_instruction(t, instruction)
        return instruction.duration
        
    def sin_ramp(self, t, duration, amplitude, offset, angular_period):
        instruction = SinRamp(t, duration, amplitude, offset, angular_period, self.parent_device.clock_limit)
        self.add_instruction(t, instruction)
        return instruction.duration
        
    def cos_ramp(self, t, duration, amplitude, offset, angular_period):
        instruction = CosRamp(t, duration, amplitude, offset, angular_period, self.parent_device.clock_limit)
        self.add_instruction(t, instruction)
        return instruction.duration
        
    def exp_ramp(self, t, duration, amplitude, offset, time_constant):
        instruction = ExpRamp(t, duration, amplitude, offset, time_constant, self.parent_device.clock_limit)
        self.add_instruction(t, instruction)
        return instruction.duration
        
class ADWinDigitalOut(DigitalOut):
    pass
    
class ADWinCard(PseudoClock):
    clock_type = 'fast clock'
    def __init__(self, name, parent_device, card_number):
        self.clock_limit = parent_device.clock_limit
        self.clock_resolution = parent_device.clock_resolution
        
        # Device must be accessed via the parent ADWin, so we must store
        # the parent's device_no as well as the card number:
        self.card_number = card_number
        self.BLACS_connection = parent_device.BLACS_connection, card_number
        # We won't call IntermediateDevice.__init__(), as we don't care
        # about the checks it does for clocking, we don't actually have
        # a clock:
        Device.__init__(self, name, parent_device, card_number)
        self.trigger_times = []
        self.wait_times = []
        self.initial_trigger_time = 0

    def trigger(self, t, *args):
        if t == 'initial':
            t = self.initial_trigger_time
            self.trigger_times.append(t)
        else:
            raise NotImplementedError("AdWins do not have waits implemented in labscript or the current firmware.")
        # Adwin cards are coordinated internally without the need for
        # triggering devices.  We split them up into pseudoclocks in
        # labscript because they are modular in nature, and it helps us
        # be compatible with AdWins that have different card setups. But
        # half of this pseudoclock stuff isn't relevant to this, so we
        # override some methods to do nothing.
        
#    def generate_code(self, hdf5_file):
#        # We don't actually need to expand out ramps and construct a pseudoclock or anything
#        # but we will anyway so that we have something to plot in runviewer
#        expanded_change_times
#        for output in self.get_all_outputs():
    
class ADWin_AO_Card(ADWinCard):
    description = 'ADWin analog output card'
    allowed_children = [AnalogOut]
    
    def generate_code(self, hdf5_file):
        Device.generate_code(self, hdf5_file)
        # This group must exist in order for BLACS to know that this
        # device is part of the experiment:
        group = hdf5_file.create_group('/devices/%s'%self.name)
        # OK, let's collect up all the analog instructions!
        self.formatted_instructions = []
        for output in self.get_all_outputs():
            for t, instruction in output.instructions.items():
                card_number = self.card_number
                channel_number = output.connection
                if isinstance(instruction, RampInstruction):
                    duration = instruction.duration
                    if isinstance(instruction, LinearRamp):
                        ramp_type = 0
                    elif isinstance(instruction, SinRamp):
                        ramp_type = 1
                    elif isinstance(instruction, CosRamp):
                        ramp_type = 2
                    elif isinstance(instruction, ExpRamp):
                        ramp_type = 3
                    else:
                        raise ValueError(instruction)
                    A = instruction.A
                    B = instruction.B
                    C = instruction.C
                else:
                    # Let's construct a ramp out of the single value instruction:
                    duration = self.clock_resolution
                    ramp_type = 0
                    A = instruction
                    B = 0
                    C = instruction
                formatted_instruction = {'t':t,
                                         'duration': duration,
                                         'card': card_number,
                                         'channel': channel_number,
                                         'ramp_type': ramp_type,
                                         'A': A, 'B': B, 'C': C}
                self.formatted_instructions.append(formatted_instruction)

                        
class ADWin_DO_Card(ADWinCard):
    description = 'ADWin digital output card'
    allowed_children = [DigitalOut]
    digital_dtype = uint32
    n_digitals = 32
    
    def generate_code(self, hdf5_file):
        Device.generate_code(self, hdf5_file)
        # This group must exist in order for BLACS to know that this
        # device is part of the experiment:
        group = hdf5_file.create_group('/devices/%s'%self.name)
        outputs = self.get_all_outputs() 
        change_times = self.collect_change_times(outputs)
        for output in outputs:
            output.make_timeseries(change_times)
        for time in change_times:
            outputarray = [0]*self.n_digitals
            for output in outputs:
                channel = output.connection
                # We have to subtract one from the channel number to get
                # the correct index, as ADWin is one-indexed, curse it.
                outputarray[channel - 1] = array(output.timeseries)
        bits = bitfield(outputarray, dtype=self.digital_dtype)
        self.formatted_instructions = []
        for t, value in zip(change_times, bits):
            formatted_instruction = {'t': t, 'card': self.card_number,'bitfield': value} 
            self.formatted_instructions.append(formatted_instruction)
            
class ADWin(PseudoClock):
    description = 'ADWin'
    allowed_children = [ADWin_AO_Card, ADWin_DO_Card]
    def __init__(self, name, device_no=1, cycle_time = default_cycle_time):
        self.BLACS_connection = device_no
        # round cycle time to the nearest multiple of 3.3333ns
        quantised_cycle_time = round(cycle_time/3.333333333333e-9)
        cycle_time = quantised_cycle_time*3.333333333333e-9
        self.clock_limit = 1./cycle_time
        self.clock_resolution = cycle_time
        Device.__init__(self, name, parent_device=None, connection=None)
        self.trigger_times = []
        self.wait_times = []
        self.initial_trigger_time = 0
    
    def do_checks(self, outputs):
        if self.trigger_times != [0]:
            raise LabscriptError('ADWin does not support retriggering or waiting.')
        for output in outputs:
            output.do_checks(self.trigger_times)

    def collect_card_instructions(self, hdf5_file):
        group = hdf5_file.create_group('/devices/%s'%self.name)
        all_analog_instructions = []
        all_digital_instructions = []
        for device in self.child_devices:
            if isinstance(device, ADWin_AO_Card):
                all_analog_instructions.extend(device.formatted_instructions)
            elif isinstance(device, ADWin_DO_Card):
                all_digital_instructions.extend(device.formatted_instructions)
            else:
                raise AssertionError("Invalid child device, shouldn't be possible")
        # Make the analog output table:
        analog_dtypes = [('t',uint), ('duration',int), ('card',int), ('channel',int),
                         ('ramp_type',int), ('A',int), ('B',int), ('C',int)]
        # sort by time:
        all_analog_instructions.sort(key=lambda instruction: instruction['t'])
        analog_data = zeros(len(all_analog_instructions)+1, dtype=analog_dtypes)
        for i, instruction in enumerate(all_analog_instructions):
            analog_data[i]['t'] = round(instruction['t']/self.clock_resolution)
            analog_data[i]['duration'] = round(instruction['duration']/self.clock_resolution)
            analog_data[i]['card'] = instruction['card']
            analog_data[i]['channel'] = instruction['channel']
            analog_data[i]['ramp_type'] = instruction['ramp_type']
            if instruction['ramp_type'] in [0]:
                # If it's a linear ramp, map the voltages for parameter A from the range [-10,10] to a uint16:
                analog_data[i]['A'] = int((instruction['A']+10)/20.*(2**16-1))
            elif instruction['ramp_type'] in [1,2,3]:
                # For an exp,  sine or cos ramp, map A from [-10,10] to a signed int16:
                analog_data[i]['A'] = int(instruction['A']/10.*(2**15-1))
            else:
                raise RuntimeError('Sanity check failed: Invalid ramp type! Something has gone wrong.')
            analog_data[i]['B'] = round(instruction['B']/self.clock_resolution) # B has units of time
            analog_data[i]['C'] = int((instruction['C']+10)/20.*(2**16-1))
        # Add the 'end of data' instruction to the end:
        analog_data[-1]['t'] = 2**32-1
        # Save to the HDF5 file:
        group.create_dataset('ANALOG_OUTS', data=analog_data)
        
        # Make the digital output table:
        digital_dtypes = [('t',uint), ('card',int), ('bitfield',int)]
        # sort by time:
        all_digital_instructions.sort(key=lambda instruction: instruction['t'])
        digital_data = zeros(len(all_digital_instructions)+1, dtype=digital_dtypes)
        for i, instruction in enumerate(all_digital_instructions):
            digital_data[i]['t'] = round(instruction['t']/self.clock_resolution)
            digital_data[i]['card'] = instruction['card']
            digital_data[i]['bitfield'] = instruction['bitfield']
        # Add the 'end of data' instruction to the end:
        digital_data[-1]['t'] = 2**32-1
        # Save to the HDF5 file:
        group.create_dataset('DIGITAL_OUTS', data=digital_data)
        group.attrs['stop_time'] = self.stop_time/self.clock_resolution
        group.attrs['cycle_time'] = self.clock_resolution
        
    def generate_code(self, hdf5_file):
        outputs = self.get_all_outputs()
        # We call the following to do the error checking it includes,
        # but we're not actually interested in the set of change times.
        # Each card will handle its own timebase issues.
        ignore = self.collect_change_times(outputs)
        self.do_checks(outputs)
        # This causes the cards to have their generate_code() methods
        # called. They collect up the instructions of their outputs,
        # and then we will collate them together into one big instruction
        # table.
        Device.generate_code(self, hdf5_file)
        self.collect_card_instructions(hdf5_file)
        
        # We don't actually care about these other things that pseudoclock
        # classes normally do, but they still do some error checking
        # that we want:
        change_times = self.collect_change_times(outputs)
        for output in outputs:
            output.make_timeseries(change_times)
        all_times, clock = self.expand_change_times(change_times, outputs)
        
