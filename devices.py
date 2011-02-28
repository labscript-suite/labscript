import control
import functions
            
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
    clock_limit = 25.0e6 # Slight underestimate I think.
    clock_connection = 11
    direct_outputs = []
    
    def perform_checks(self):
        for output in self.outputs:
            if output.connected_to_device is self:
                if not (isinstance(output.connection_number,int) or output.connection_number < 12):
                    sys.stderr.write('%s is set as connected to output connection %d of %s. Output connection \
                                      number must be a integer less than 12\n.'%(output.name, output.connection_number, self.name))
                    sys.exit(1)
                for other_output in self.direct_outputs:
                    if output.connection_number == other_output.connection_number:
                        sys.stderr.write('%s %s and %s %s are both set as connected to output %d of %s! Stopping.\n'%(output.name,
                                         other_output.name, output.connection_number, self.name))
                        sys.exit(1)
                self.direct_outputs.append(output)    
                
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
        # Gotta put a stop instruction at the end. It will have a short
        # delay time and set everything back to zero:
        self.pb_inst.append({'flags': '000000000000', 'instruction': 'STOP',
                                 'data': 0, 'delay': 10.0/self.clock_limit*1e9})  
                                           
    def write_instructions_to_files(self):
        # The raw output of each analogue output not direcly connected
        # to the PulseBlaster:
        for output in self.outputs:
            if output not in self.direct_outputs:
                output.write_raw_output_to_file()
        # The table of instructions for the PulseBlaster itself:
        with open(self.name+'.dat','w') as outfile:
            for inst in self.pb_inst:
                flagint = '%04d'%int(inst['flags'][::-1],2)
                instructionint = str(self.pb_instructions[inst['instruction']])
                dataint = '%04d'%inst['data']
                delaydouble = repr(inst['delay']) # repr to keep high precision
                outfile.write('\t'.join(['0']*10 + [flagint,instructionint,dataint,delaydouble,'\n']))
        print 'saved', self.name+'.dat'
             
    def generate_code(self):
        self.perform_checks()
        self.make_instruction_table()
        self.convert_to_pb_inst()
        self.write_instructions_to_files()


class AnalogueOut(control.Output):
    description = 'analogue output'
    def ramp(self,t,duration,initial,final,samplerate):
        self.add_instruction(t, {'function': functions.ramp(t,duration,initial,final), 'description':'linear ramp',
                                 'end time': t + duration, 'clock rate': samplerate})
                                 
    def sine(self,t,duration,amplitude,angfreq,phase,dc_offset,samplerate):
        self.add_instruction(t, {'function': functions.sine(t,amplitude,angfreq,phase,dc_offset), 'description':'sine wave',
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
         
                               

    
