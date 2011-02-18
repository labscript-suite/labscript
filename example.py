from control import PulseBlaster, Shutter, plot_timeseries, run, stop
import time
import profile

pulseblaster_one = PulseBlaster('pulseblaster_one')
shutter_one = Shutter('shutter_one',pulseblaster_one,0)
shutter_two = Shutter('shutter_two',pulseblaster_one,1)

shutter_one.close(0)
shutter_two.close(0)
#shutter_two.open(0)

for t in range(1,10):
    shutter_one.open(t)
    shutter_one.close((t+0.5))
for t in range(10,20):
    shutter_two.open(t)
    shutter_two.close((t+0.5))

run()
time.sleep(20)
stop()
plot_timeseries()

