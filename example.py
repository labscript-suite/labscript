import time
import shutil
start_time = time.time()
very_start_time = time.time()
import os
print "other imports", time.time() - start_time
start_time = time.time()
shutil.copy('dummy_template.h5', 'dummy.h5')
print "copying the dummy file:", time.time() - start_time

start_time = time.time()
from labscript import *
print "from labscript import *", time.time() - start_time
start_time = time.time()
print z, A, B, C, D

pulseblaster1 = PulseBlaster('PulseBlaster',stop_time=11)
NI_board1 = NIBoard('NI PCI-6733', pulseblaster1)
devices = []

analogue1 = AnalogueOut('output 1', NI_board1,'AO0')
analogue2 = AnalogueOut('output 2', NI_board1,'AO1')
analogue3 = AnalogueOut('output 3', NI_board1,'AO2')

shutter1 = Shutter('shutter 1', NI_board1, 'D01')
shutter1.close(t=0)
shutter1.open(t=5.89)
analogue1.constant(t=0,value=2)
analogue1.ramp(t=1, duration=2, initial=2, final=3, samplerate=0.5e6)

analogue2.constant(t=0,value=3)
analogue2.ramp(t=2, duration=3, initial=3, final=4, samplerate=0.5e6)
analogue2.constant(5.9,5)
analogue2.constant(7,4)
analogue2.constant(8,5)
analogue3.sine(t=0,duration=10,amplitude=10,angfreq=2,phase=0,dc_offset=0.0,samplerate=0.5e6)
print "actual script", time.time() - start_time
start_time = time.time()
generate_code()
print "pulseblaster1.generate_code()", time.time() - start_time
plot_outputs(display=True)
start_time = time.time()
os.system('sync') # linux only to measure hard drive write time, which is otherwise deferred.
print "os.system(sync)", time.time() - start_time
print "total time:", time.time() - very_start_time

