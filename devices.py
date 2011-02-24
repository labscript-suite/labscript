import control
import time

class NIBoard:
    description = 'NI board'
    clock_limit = 1e6 #TODO get a real number to put here
    def __init__(self, name, clockedby):
        self.name = name
        # Make sure the clock knows is limited by this device's max clock rate:
        if clockedby.clock_limit > self.clock_limit:
            clockedby.clock_limit = self.clock_limit
        self.outputs = clockedby.outputs
       
class PulseBlaster(control.IODevice):
    description = 'PulseBlaster'
    clock_limit = 25.0e6 # TODO get a real number to put here
    clock_connection = 11
    direct_outputs = []
    def generate_code(self):
        self.make_instruction_table()
        for output in self.outputs:
            if output.connected_to_device is self:
                self.direct_outputs.append(output)
        print 'direct outputs are:', [out.name for out in self.direct_outputs]
 
class AnalogueOut(control.Output):
    description = 'analogue output'
    def ramp(self,t,duration,initial,final,samplerate):
        self.add_instruction(t, {'function': control.ramp(t,duration,initial,final), 
                                 'end time': t + duration, 'clock rate': samplerate})
                                 
    def sine(self,t,duration,amplitude,angfreq,phase,dc_offset,samplerate):
        self.add_instruction(t, {'function': control.sine(t,amplitude,angfreq,phase,dc_offset), 
                                 'end time': t + duration, 'clock rate': samplerate})
                                 
    def constant(self,t,value):
        self.add_instruction(t,value)
        
class DigitalOut(control.Output):
    description = 'digital output'
    allowed_states = {1:'high', 0:'low'}
    def go_high(self,t):
        self.add_instruction(t,1)
    def go_low(self,t):
        self.add_instruction(t,0) 

class Shutter(DigitalOut):
    description = 'shutter'
    allowed_states = {1:'open', 0:'closed'}  
    def open(self,t):
        self.go_high(t)
    def close(self,t):
        self.go_low(t)
                               
if __name__ == '__main__':
    pulseblaster1 = PulseBlaster('PulseBlaster_1')
    NI_board1 = NIBoard('NI_board_1', pulseblaster1)
    analogue1 = AnalogueOut('output 1', NI_board1,0)
    analogue2 = AnalogueOut('output 2', NI_board1,1)
    analogue3 = AnalogueOut('output 3', NI_board1,2)
    
    shutter1 = Shutter('shutter 1', pulseblaster1,0)
    shutter1.close(t=0)
    shutter1.open(t=5.89)
    analogue1.constant(t=0,value=2)
    analogue1.ramp(t=1, duration=2, initial=2, final=3, samplerate=5)

    analogue2.constant(t=0,value=3)
    analogue2.ramp(t=2, duration=3, initial=3, final=4, samplerate=10)
    analogue2.constant(5.9,5)
    analogue2.constant(7,4)
    analogue2.constant(8,5)
    analogue3.sine(t=0,duration=10,amplitude=1,angfreq=2,phase=0,dc_offset=0.5,samplerate=3)

    pulseblaster1.generate_code()
    control.plot_outputs()
    
