import shutil
shutil.copy('dummy_template.h5', 'dummy.h5')

from labscript import *

pulseblaster1 = PulseBlaster(1)
NI_board1 = NI_PCIe_6363(1, pulseblaster1,'fast clock',acquisition_rate=1e3)
novatech1 = NovaTechDDS9M(1, pulseblaster1,'slow clock')

analog0 = AnalogOut('output 1', NI_board1,'ao0')
analog1 = AnalogOut('output 2', NI_board1,'ao1')
analog2 = AnalogOut('output 3', NI_board1,'ao2')
input1 = AnalogIn('input 1', NI_board1,'ai0')

shutter1 = Shutter('shutter 1', NI_board1, 'port0/line0', delay='calibrated')
shutter2 = Shutter('shutter 2', pulseblaster1, 'flag 2', delay='calibrated')
dds1 = DDS('DDS 1', novatech1, 'channel 0')
dds2 = DDS('DDS 2', novatech1, 'channel 1')

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
