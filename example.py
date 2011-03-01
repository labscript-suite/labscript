import os
os.system('python dummy.py')

from labscript import *

print 'params are:'
for key,val in params.items():
    print str(key) + ': ' + str(val)

pulseblaster1 = PulseBlaster('PulseBlaster',stop_time=11)
NI_board1 = NIBoard('NI PCI-6733', pulseblaster1)

analogue1 = AnalogueOut('output 1', NI_board1,'AO0')
analogue2 = AnalogueOut('output 2', NI_board1,'AO1')
analogue3 = AnalogueOut('output 3', NI_board1,'AO2')

shutter1 = Shutter('flag 1', pulseblaster1,0)
shutter1.close(t=0)
shutter1.open(t=5.89)
analogue1.constant(t=0,value=2)
analogue1.ramp(t=1, duration=2, initial=2, final=3, samplerate=0.5e6)

analogue2.constant(t=0,value=3)
analogue2.ramp(t=2, duration=3, initial=3, final=4, samplerate=0.5e6)
analogue2.constant(5.9,5)
analogue2.constant(7,4)
analogue2.constant(8,5)
analogue3.sine(t=0,duration=10,amplitude=1,angfreq=2,phase=0,dc_offset=0.5,samplerate=0.5e6)

pulseblaster1.generate_code()
#plot_outputs()
