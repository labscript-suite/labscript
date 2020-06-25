#####################################################################
#                                                                   #
# /example.py                                                       #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

import numpy as np
from labscript import (
    AnalogIn,
    AnalogOut,
    ClockLine,
    DDS,
    DigitalOut,
    MHz,
    Shutter,
    StaticDDS,
    WaitMonitor,
    start,
    stop,
    wait,
)
from labscript_devices.PulseBlaster import PulseBlaster
from labscript_devices.NI_DAQmx.labscript_devices import NI_PCIe_6363
from labscript_devices.NovaTechDDS9M import NovaTechDDS9M
from labscript_devices.Camera import Camera
from labscript_devices.PineBlaster import PineBlaster
from labscript_devices.NI_DAQmx.labscript_devices import NI_PCI_6733
from labscript_utils.unitconversions import example1, example2, example3

PulseBlaster(name="pulseblaster_0", board_number=0)
ClockLine(
    name="pulseblaster_0_clockline_fast",
    pseudoclock=pulseblaster_0.pseudoclock,
    connection="flag 0",
)
ClockLine(
    name="pulseblaster_0_clockline_slow",
    pseudoclock=pulseblaster_0.pseudoclock,
    connection="flag 1",
)
NI_PCIe_6363(
    name="ni_card_0",
    parent_device=pulseblaster_0_clockline_fast,
    clock_terminal="ni_pcie_6363_0/PFI0",
    MAX_name="ni_pcie_6363_0",
    acquisition_rate=100e3,
)
NovaTechDDS9M(
    name="novatechdds9m_0",
    parent_device=pulseblaster_0_clockline_slow,
    com_port="com10",
)

# Create a BIAS Camera, tirggered to take photos with flag 3 of pulseblaster_0
Camera(
    name="andor_ixon_0",
    parent_device=pulseblaster_0.direct_outputs,
    connection="flag 3",
    BIAS_port=42520,
    serial_number="0000",
    SDK="IMAQdx",
    effective_pixel_size=4.6e-6,
    exposure_time=0.1,
    orientation="top",
)

# A second pseudoclock to just clock a NI_PCI_6733 Card
PineBlaster(
    name="pineblaster_0",
    trigger_device=ni_card_0,
    trigger_connection="port0/line15",
    usbport="COM7",
)
NI_PCI_6733(
    name="ni_card_1",
    parent_device=pineblaster_0.clockline,
    clock_terminal="ni_pcie_6733_0/PFI0",
    MAX_name="ni_pci_6733_0",
)

# Create the output/input channels on the above devices use the example1 conversion
# class located in pythonlib/unitconversions/example.py with default paremeters
AnalogOut(
    name="analog0",
    parent_device=ni_card_1,
    connection="ao0",
    unit_conversion_class=example1,
)

# same as above, but we are changing some parameters used in the conversion and
# specifying a prefix to be used with units. You can now program in mA, uA, mGauss,
# uGauss
AnalogOut(
    name="analog1",
    parent_device=ni_card_1,
    connection="ao1",
    unit_conversion_class=example1,
    unit_conversion_parameters={"a": 5, "b": 1, "magnitudes": ["m", "u"]},
)
AnalogOut(name="analog2", parent_device=ni_card_0, connection="ao2")
AnalogOut(name="analog3", parent_device=ni_card_0, connection="ao3")
AnalogIn(name="input1", parent_device=ni_card_0, connection="ai0")
DigitalOut(name="do0", parent_device=ni_card_0, connection="port0/line2")
Shutter(
    name="shutter1", parent_device=ni_card_0, connection="port0/line1", delay=(0, 0)
)
Shutter(
    name="shutter2",
    parent_device=pulseblaster_0.direct_outputs,
    connection="flag 2",
    delay=(0, 0),
)
DigitalOut(
    name="switch", parent_device=pulseblaster_0.direct_outputs, connection="flag 4"
)

DDS(name="dds1", parent_device=novatechdds9m_0, connection="channel 0")
DDS(name="dds2", parent_device=novatechdds9m_0, connection="channel 1")
StaticDDS(name="dds5", parent_device=novatechdds9m_0, connection="channel 2")
# The next DDS is special because it has the frequency and amplitude calibrated using
# example2 and example3 classes from pythonlib/unitconversions/example.py
DDS(
    name="dds3",
    parent_device=pulseblaster_0.direct_outputs,
    connection="dds 0",
    freq_conv_class=example2,
    freq_conv_params={"a": 4, "b": 6},
    amp_conv_class=example3,
    amp_conv_params={"a": 2, "b": 22, "magnitudes": ["m"]},
)
DDS(name="dds4", parent_device=pulseblaster_0.direct_outputs, connection="dds 1")

# This sets up the inputs/counters/etc that will monitor
# The first paremeter is the name for the WaitMonitor device
# The second and third paremeters are the device and channel respectively that goes
# high when a wait begins and low when it ends. This output should be
# physically connected to a counter specified in the next two paremeters.
# The final two paremeters specify the device/channel that is to trigger the
# pseudoclock if the WAIT instruction hits the specified timeout. The output of
# this channel should be physicaly connect to the external trigger of the master
# pseudoclock.
WaitMonitor(
    name="wait_monitor",
    parent_device=ni_card_0,
    connection="port0/line0",
    acquisition_device=ni_card_0,
    acquisition_connection="ctr0",
    timeout_device=ni_card_0,
    timeout_connection="pfi1",
)

# A variable to define the acquisition rate used for the analog outputs below.
# This is just here to show you that you can use variables instead of typing in numbers!
# Furthermore, these variables could be defined within runmanager (rather than in the
# code like this one is)
# for easy manipulation via a graphical interface.
rate = 1e4

# The time (in seconds) we wish the pineblaster pseudoclock to start after
# the master pseudoclock (the pulseblaster)
pineblaster_0.set_initial_trigger_time(t=0.9)

# Start the experiment!
start()

# A variable to keep track of time
t = 0

# Analog Acquisitions are acquired at the sample rate specified when the *device* is
# instantiated (eg NI_PCIE_6363() above)
# Acquire an analog trace on this channel from t=0s to t=1s
input1.acquire(label='measurement1', start_time=0, end_time=1)
# Acquire an analog trace on this channel from t=3s to t=5s
input1.acquire(label='measurement2', start_time=3, end_time=5)
# Acquire an analog trace on this channel from t=7s to t=9s
input1.acquire(label='measurement3', start_time=7, end_time=9)

# Set some initial values (t=0) for DDS 1
dds1.setamp(t=t, value=0.5)
dds1.setfreq(t=t, value=0.6)
dds1.setphase(t=t, value=0.7)

# Set some values for dds2 at t=1s. They will have value '0' before this
# time unless otherwise set
dds2.setamp(t=t + 1, value=0.9)
dds2.setfreq(t=t + 1, value=1.0)
dds2.setphase(t=t + 1, value=1.1)

# dds5 is a "static" DDS. This means its value can only be set once, and
# will be set just before the experiment begins
dds5.setfreq(value=90 * MHz)
dds5.setamp(value=1)
dds5.setphase(value=0)

# Have the shutters start in the closed state (t=0)
shutter1.close(t=t)
shutter2.close(t=t)

# Analog0 is attached to ni_card_1, which is an NI-pci_6733 card (MAX name
# ni_pcie_6733_0) clocked by a pineblaster
# The pineblaster is being triggered to start by a pulseblaster, which introduces
# a delay into the start of output from the NI card
# This is all handled by labscript, and you still specify times from the beginning of
# the experiment (when the master pseudoclock is started)
# YOU DO NOT HAVE TO TAKE INTO ACCOUNT THE DELAY YOURSELF!!
# You do however need to make sure you do not command output from this device before
# the device has actually started.
# To do so, make sure no commands happen on this channel before analog0.t0
# (this variable contains the delay time!)
analog0.constant(t=analog0.t0, value=2)

# Set this channel to a constant value
# As this channel is clocked by the master pseudoclock, you can command
# output from t=0
analog2.constant(t=t, value=3)

# Again, this must not start until analog1.t0 or later!
analog1.sine(
    t=analog1.t0,
    duration=6,
    amplitude=5,
    angfreq=2 * np.pi,
    phase=0,
    dc_offset=0.0,
    samplerate=rate,
)

# Let's increment our time variable!
t += max(1, analog0.t0)

# Open the shutter, enable the DDS, ramp an analog output!
shutter2.open(t=t)
dds3.enable(t=t)
analog0.ramp(t=t, duration=2, initial=2, final=3, samplerate=rate)

# Take a picture
andor_ixon_0.expose(name='exposure_1', t=t, frametype='flat')
andor_ixon_0.expose(name='exposure_1', t=t + 1, frametype='atoms')

# Do some more things at various times!
# (these are ignoring the t variable)
def my_ramp(t, *args, **kwargs):
    lambda_func = functions.sine_ramp(*args, **kwargs)
    return lambda_func(t)


analog2.sine_ramp(
    t=2.25, duration=3, initial=3, final=4, samplerate=rate, truncation=0.7
)
shutter1.open(t=5.89)
analog2.constant(t=5.9, value=5)
analog2.constant(t=7, value=4)
analog2.constant(t=8, value=5)

# Incremenent t by 9 seconds
t += 9

# Wait for an external trigger on the master pseudoclock
# Waits must be names
# The timeout defaults to 5s, unless otherwise specified.
# The timeout specifies how long to wait without seeing the external
# trigger before continuing the experiment
t += wait(label='my_first_wait', t=t, timeout=2)

# Waits take very little time as far as labscript is concerned. They only add on the
# retirggering time needed to start devices up and get them all in sync again.
# After a wait, labscript time (the t variable here) and execution time (when the
# hardware instructions are executed on the hardware) will not be the same
# as the wait instruction may take anywhere from 0 to "timeout" seconds,
# and this number is only determined during execution time.

t += 1
# Do something 1s after the wait!
switch.go_high(t=t)

# Examples programming in different units as specified in the
# unitconversion classes passed as parameters to the channel definition
analog0.constant(t=t, value=5, units='A')
analog1.constant(t=t, value=1000, units='mGauss')
dds3.setfreq(t=t, value=50, units='detuned_MHz')
dds3.setamp(t=t, value=1.9, units='W')

# Hold values for 2 seconds
t += 2

analog0.ramp(t=t, duration=1, initial=5, final=7, samplerate=rate, units='Gauss')
analog1.constant(t=t, value=3e6, units='uA')
dds3.setfreq(t=t, value=60, units='detuned_MHz')
dds3.setamp(t=t, value=500, units='mW')

# Hold values for 2 seconds
t += 2

# Stop at t=15 seconds, note that because of the wait timeout, this might
# be as late as 17s (Plus a little bit of retriggering time) in execution
# time
stop(t=t)
