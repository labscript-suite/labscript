from labscript import *

PulseBlaster(  'pulseblaster_0')
NI_PCIe_6363(  'ni_pcie_6363_0',  pulseblaster_0, 'fast clock', 'ni_pcie_6363_0/PFI0')
NI_PCI_6733(   'ni_pci_6733_0',  pulseblaster_0, 'fast clock', 'ni_pcie_6363_0/PFI0')
NovaTechDDS9M( 'novatechdds9m_0', pulseblaster_0, 'slow clock')

AnalogOut( 'analog0',  ni_pcie_6363_0,         'ao0')
AnalogOut( 'analog1',  ni_pcie_6363_0,         'ao1')
AnalogOut( 'analog2',  ni_pcie_6363_0,         'ao2')
AnalogIn(   'input1',  ni_pcie_6363_0,         'ai0')
Shutter(  'shutter1',  ni_pcie_6363_0, 'port0/line0', delay=(0,0))
Shutter(  'shutter2',  pulseblaster_0,      'flag 2', delay=(0,0))
DDS(          'dds1', novatechdds9m_0,   'channel 0')
DDS(          'dds2', novatechdds9m_0,   'channel 1')
StaticDDS(    'dds5', novatechdds9m_0,   'channel 2')
DDS(          'dds3',  pulseblaster_0,       'dds 0')
DDS(          'dds4',  pulseblaster_0,       'dds 1')
Camera('andor_ixon_0', pulseblaster_0,   'flag 3', exposuretime=.1,orientation='top')

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

dds5.setamp(1)
shutter1.close(t)
shutter2.close(t)
analog0.constant(t,2)

analog2.constant(t,3)
analog1.sine(t,duration=6*scale,amplitude=5,angfreq=2*pi,phase=0,dc_offset=0.0,samplerate=rate)
t = 1*scale
shutter2.open(t)
dds3.enable(t)
analog0.ramp(t, duration=2*scale, initial=2, final=3, samplerate=rate)

andor_ixon_0.expose('exposure_1',t,'flat')
andor_ixon_0.expose('exposure_1',t+1,'atoms')

analog2.ramp(t=2*scale, duration=3*scale, initial=3, final=4, samplerate=rate)
shutter1.open(t=5.89*scale)
analog2.constant(t=5.9*scale,value=5)
analog2.constant(t=7*scale,value=4)
analog2.constant(t=8*scale,value=5)

analysis(['av_measurement_1'], 'pylab', 'average', traces=['measurement_1'], globals_=[], other_results = [], args = [], kwargs = {})
analysis(['sqrt_av_measurement_1'], 'math','sqrt', other_results=['av_measurement_1'] )

stop(t=8*scale+2e-6)

















