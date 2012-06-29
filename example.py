from labscript import *
from unitconversions import *

RFBlaster('rfblaster_0','some.ip.address')
DDS(          'dds3',  rfblaster_0,       'dds 0',freq_conv_class=test,freq_conv_params={'a':4,'b':6},amp_conv_class=test,amp_conv_params={'a':2,'b':22})
DDS(          'dds4',  rfblaster_0,       'dds 1')

t = 0

dds3.setamp(t,1)
dds3.setfreq(t,1)
dds4.frequency.ramp(t,duration=10,initial=0,final=5,samplerate=10)
stop(11)









