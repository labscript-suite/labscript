import shutil, os
shutil.copy(os.path.join('utils','dummy_template.h5'), 'dummy.h5')

from labscript import *

PulseBlaster(  'pulseblaster_0')
NI_PCIe_6363(  'ni_pcie_6363_0',  pulseblaster_0, 'fast clock', acquisition_rate=1e3)
NovaTechDDS9M( 'novatechdds9m_0', pulseblaster_0, 'slow clock')

AnalogOut( 'analog0',  ni_pcie_6363_0,         'ao0')
AnalogOut( 'analog1',  ni_pcie_6363_0,         'ao1')
AnalogOut( 'analog2',  ni_pcie_6363_0,         'ao2')
AnalogIn(   'input1',  ni_pcie_6363_0,         'ai0')
Shutter(  'shutter1',  ni_pcie_6363_0, 'port0/line0', delay='calibrated')
Shutter(  'shutter2',  pulseblaster_0,      'flag 2', delay='calibrated')
DDS(          'dds1', novatechdds9m_0,   'channel 0')
DDS(          'dds2', novatechdds9m_0,   'channel 1')

scale = 1.0
rate = 1e4
t = 0

input1.acquire('measurement1',0*scale,1*scale)
input1.acquire('measurement2',3*scale,5*scale)
input1.acquire('measurement3',7*scale,9*scale)

dds1.setamp(t,0.5)
dds1.setfreq(t,0.6)
dds1.setphase(t,0.7)
dds2.setamp(t,0.8)
dds2.setfreq(t,0.9)
dds2.setphase(t,1.0)

shutter1.close(t)
shutter2.close(t)
analog0.constant(t,2)
analog2.constant(t,3)
analog1.sine(t,duration=10*scale,amplitude=5,angfreq=2*pi,phase=0,dc_offset=0.0,samplerate=rate)
t = 1*scale
shutter2.open(t)
analog0.ramp(t, duration=2*scale, initial=2, final=3, samplerate=rate)

analog2.ramp(t=2*scale, duration=3*scale, initial=3, final=4, samplerate=rate)
shutter1.open(t=5.89*scale)
analog2.constant(t=5.9*scale,value=5)
analog2.constant(t=7*scale,value=4)
analog2.constant(t=8*scale,value=5)

stop(t=10*scale)

plot_outputs()
