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
    pb_instructions = {'STOP': 1, 'LOOP': 2, 'END_LOOP': 3}
    description = 'PulseBlaster'
    clock_limit = 25.0e6 # TODO get a real number to put here
    clock_connection = 11
    direct_outputs = []
    
    def convert_to_pb_inst(self):
        self.pb_inst = []
        # index to keep track of where in output.raw_output the
        # pulseblaster flags are coming from
        i = 0
        # index to record what line number of the pulseblaster hardware
        # instructions we're up to:
        j = 0
        for instruction in self.clock:
            flags = [0]*12
            for output in self.direct_outputs:
                flags[output.connection_number] = int(output.raw_output[i])
            flags[11] = 1
            flagstring = ''.join([str(flag) for flag in flags])
            if instruction['reps'] > 1048576:
                raise Exception('Pulseblaster cannot support more than 1048576 loop iterations. ' +
                                 str(instruction['reps']) +' were requested at t = ' + str(instruction['start']) + '. '+
                                 'This can be fixed easily enough by using nested loops. If it is needed, ' +
                                  'please file a feature request at http://redmine.physics.monash.edu.au/projects/labscript.')
            self.pb_inst.append({'flags': flagstring, 'instruction': 'LOOP',
                                 'data': instruction['reps'], 'delay': instruction['step']*1e9/2.0})
            flags[11] = 0
            flagstring = ''.join([str(flag) for flag in flags])
            self.pb_inst.append({'flags': flagstring, 'instruction': 'END_LOOP',
                                 'data': j, 'delay': instruction['step']*1e9/2.0})
            j += 2
            i += instruction['reps']
        for thing in self.pb_inst:
            print thing['flags'],thing['instruction'],thing['data'], thing['delay']
            
    def generate_code(self):
        self.make_instruction_table()
        for output in self.outputs:
            if output.connected_to_device is self:
                if not (isinstance(output.connection_number,int) or output.connection_number < 12):
                    raise ValueError('%s is set as connected to output connection %d of %s. Output connection \
                                      number must be a integer less than 12.'%(output.name, output.connection_number, self.name))
                for other_output in self.direct_outputs:
                    if output.connection_number == other_output.connection_number:
                        raise ValueError('%s %s and %s %s are both set as connected to output %d of %s!'%(output.name,
                                         other_output.name, output.connection_number, self.name))
                self.direct_outputs.append(output)
        self.convert_to_pb_inst()
                
#        print 'direct outputs are:\n', '\n'.join([out.name + ' connected to output ' + str(out.connection_number) for out in self.direct_outputs])
#        print 'start','   step','    reps'
#        ticks = 0
#        for thing in self.clock:
#            ticks += thing['reps']
#            if isinstance(thing, dict):
#                print '%1f'%thing['start'], '%1f'%thing['step'], '%02d'%thing['reps']
#            else:
#                print round(thing,2)
#        print 'total no of clock ticks:', ticks
#        print 'total no of timepoints:', len(self.flat_times)


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
    time.sleep(1)
    pulseblaster1 = PulseBlaster('PulseBlaster_1',stop_time=11)
    NI_board1 = NIBoard('NI_board_1', pulseblaster1)
    analogue1 = AnalogueOut('output 1', NI_board1,0)
    analogue2 = AnalogueOut('output 2', NI_board1,1)
    analogue3 = AnalogueOut('output 3', NI_board1,2)
    
    #shutter1 = Shutter('shutter 1', pulseblaster1,0)
    #shutter1.close(t=0)
    #shutter1.open(t=5.89)
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
    
