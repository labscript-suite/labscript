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

import __init__ # only have to do this because we're inside the labscript directory
from labscript import *
from labscript_devices.PulseBlaster import PulseBlaster
from labscript_devices.NI_PCIe_6363 import NI_PCIe_6363
from labscript_devices.NovaTechDDS9M import NovaTechDDS9M
from labscript_devices.Camera import Camera
from labscript_devices.PineBlaster import PineBlaster
from labscript_devices.NI_PCI_6733 import NI_PCI_6733
from labscript_utils.unitconversions import *

PulseBlaster(name='pulseblaster_0', board_number=0)
ClockLine(name='pulseblaster_0_clockline_fast', pseudoclock=pulseblaster_0.pseudoclock, connection='flag 0')
ClockLine(name='pulseblaster_0_clockline_slow', pseudoclock=pulseblaster_0.pseudoclock, connection='flag 1')
NI_PCIe_6363(name='ni_card_0',  parent_device=pulseblaster_0_clockline_fast, clock_terminal='ni_pcie_6363_0/PFI0', MAX_name='ni_pcie_6363_0', acquisition_rate=100e3)
NovaTechDDS9M(name='novatechdds9m_0', parent_device=pulseblaster_0_clockline_slow, com_port="com10")

# Create a BIAS Camera, tirggered to take photos with flag 3 of pulseblaster_0
Camera('andor_ixon_0', pulseblaster_0.direct_outputs,   'flag 3',BIAS_port = 42520,serial_number="0000", SDK="IMAQdx", effective_pixel_size = 4.6e-6, exposuretime=.1, orientation='top')

# A second pseudoclock to just clock a NI_PCI_6733 Card
PineBlaster(name='pineblaster_0', trigger_device=ni_card_0, trigger_connection='port0/line15', usbport='COM7')
NI_PCI_6733(name='ni_card_1',  parent_device=pineblaster_0.clockline, clock_terminal='ni_pcie_6733_0/PFI0', MAX_name='ni_pci_6733_0')

# Create the output/input channels on the above devices
AnalogOut( 'analog0',  ni_card_1,         'ao0',         unit_conversion_class=example1) # use the example1 conversion class located in pythonlib/unitconversions/example.py with default paremeters

# same as above, but we are changing some parameters used in the conversion and specifying a prefix to be used with units. You can now program in mA, uA, mGauss, uGauss
AnalogOut( 'analog1',  ni_card_1,         'ao1',         unit_conversion_class=example1,unit_conversion_parameters = {'a':5,'b':1, 'magnitudes':['m','u']}) 
AnalogOut( 'analog2',  ni_card_0,         'ao2')
AnalogIn(   'input1',  ni_card_0,         'ai0')
Shutter(  'shutter1',  ni_card_0,         'port0/line1', delay=(0,0))
Shutter(  'shutter2',  pulseblaster_0.direct_outputs,    'flag 2',      delay=(0,0))
DigitalOut( 'switch',  pulseblaster_0.direct_outputs,    'flag 4')

DDS(          'dds1',  novatechdds9m_0,   'channel 0')
DDS(          'dds2',  novatechdds9m_0,   'channel 1')
StaticDDS(    'dds5',  novatechdds9m_0,   'channel 2')
# The next DDS is special because it has the frequency and amplitude calibrated using example2 and example3 classes from pythonlib/unitconversions/example.py
DDS(          'dds3',  pulseblaster_0.direct_outputs,    'dds 0',       freq_conv_class=example2,freq_conv_params={'a':4,'b':6},amp_conv_class=example3,amp_conv_params={'a':2,'b':22,'magnitudes':['m']})
DDS(          'dds4',  pulseblaster_0.direct_outputs,    'dds 1')

# This sets up the inputs/counters/etc that will monitor
# The first paremeter is the name for the WaitMonitor device
# The second and third paremeters are the device and channel respectively that goes high when a wait begins and low when it ends. This output should be   
# physically connected to a counter specified in the next two paremeters.
# The final two paremeters specify the device/channel that is to trigger the pseudoclock if the WAIT instruction hits the specified timeout. The output of 
# this channel should be physicaly connect to the external trigger of the master pseudoclock.
WaitMonitor('wait_monitor', ni_card_0, 'port0/line0', ni_card_0, 'ctr0', ni_card_0, 'pfi1')

# A variable to define the acquisition rate used for the analog outputs below.
# This is just here to show you that you can use variables instead of typing in numbers!
# Furthermore, these variables could be defined within runmanager (rather than in the code like this one is)
# for easy manipulation via a graphical interface.
rate = 1e4

# The time (in seconds) we wish the pineblaster pseudoclock to start after the master pseudoclock (the pulseblaster)
pineblaster_0.set_initial_trigger_time(1)

# Start the experiment!
start()

# A variable to keep track of time
t = 0

# Analog Acquisitions are acquired at the sample rate specified when the *device* is instantiated (eg NI_PCIE_6363() above)
# Acquire an analog trace on this channel from t=0s to t=1s
input1.acquire('measurement1',0,1)
# Acquire an analog trace on this channel from t=3s to t=5s
input1.acquire('measurement2',3,5)
# Acquire an analog trace on this channel from t=7s to t=9s
input1.acquire('measurement3',7,9)

# Set some initial values (t=0) for DDS 1
dds1.setamp(t,0.5)
dds1.setfreq(t,0.6)
dds1.setphase(t,0.7)

# Set some values for dds2 at t=1s. They will have value '0' before this time unless otherwise set
dds2.setamp(t+1,0.9)
dds2.setfreq(t+1,1.0)
dds2.setphase(t+1,1.1)

# dds5 is a "static" DDS. This means its value can only be set once, and will be set just before the experiment begins
dds5.setfreq(90*MHz)
dds5.setamp(1)
dds5.setphase(0)

# Have the shutters start in the closed state (t=0)
shutter1.close(t)
shutter2.close(t)

# Analog0 is attached to ni_card_1, which is an NI-pci_6733 card (MAX name ni_pcie_6733_0) clocked by a pineblaster
# The pineblaster is being triggered to start by a pulseblaster, which introduces a delay into the start of output from the NI card
# This is all handled by labscript, and you still specify times from the beginning of the experiment (when the master pseudoclock is started)
# YOU DO NOT HAVE TO TAKE INTO ACCOUNT THE DELAY YOURSELF!!
# You do however need to make sure you do not command output from this device before the device has actually started.
# To do so, make sure no commands happen on this channel before analog0.t0 (this variable contains the delay time!)
analog0.constant(analog0.t0,2)

# Set this channel to a constant value
# As this channel is clocked by the master pseudoclock, you can command output from t=0
analog2.constant(t,3)

# Again, this must not start until analog1.t0 or later!
analog1.sine(analog1.t0,duration=6,amplitude=5,angfreq=2*pi,phase=0,dc_offset=0.0,samplerate=rate)

# Let's increment our time variable!
t += max(1,analog0.t0)

# Open the shutter, enable the DDS, ramp an analog output!
shutter2.open(t)
dds3.enable(t)
analog0.ramp(t, duration=2, initial=2, final=3, samplerate=rate)

# Take a picture
andor_ixon_0.expose('exposure_1',t,'flat')
andor_ixon_0.expose('exposure_1',t+1,'atoms')

# Do some more things at various times!
# (these are ignoring the t variable)
analog2.ramp(t=2.25, duration=3, initial=3, final=4, samplerate=rate)
shutter1.open(t=5.89)
analog2.constant(t=5.9,value=5)
analog2.constant(t=7,value=4)
analog2.constant(t=8,value=5)


t+=9 # set t=10 seconds
# Wait for an external trigger on the master pseudoclock
# Waits must be names
# The timeout defaults to 5s, unless otherwise specified.
# The timeout specifies how long to wait without seeing the external trigger before continuing the experiment 
t += wait('my_first_wait', t=t, timeout=2)
# Waits take very little time as far as labscript is concerned. They only add on the retirggering time needed to start devices up and get them all in sync again.
# After a wait, labscript time (the t variable here) and execution time (when the hardware instructions are executed on the hardware) will not be the same 
# as the wait instruction may take anywhere from 0 to "timeout" seconds, and this number is only determined during execution time.

t+=1
# Do something 1s after the wait!
switch.go_high(t)

# Examples programming in different units as specified in the unitconversion classes passed as parameters to the channel definition
analog0.constant(t, value=5, units='A')
analog1.constant(t, value=1000, units='mGauss')
dds3.setfreq(t,value=50,units='detuned_MHz')
dds3.setamp(t,value=1.9,units='W')
t+=2
analog0.ramp(t,duration=1,initial=5,final=7,samplerate=rate, units='Gauss')
analog1.constant(t, value=3e6, units='uA')
dds3.setfreq(t,value=60,units='detuned_MHz')
dds3.setamp(t,value=500,units='mW')
t+=2
# Stop at t=15 seconds, note that because of the wait timeout, this might be as late as 17s (Plus a little bit of retriggering time) in execution time
stop(t)
