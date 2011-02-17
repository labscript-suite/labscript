
class IODevice:
    type_string = 'generic input/output device'
    n_outputs = {'analogue':None,'digital':None}
    n_inputs = {'analogue':None,'digital':None}
    
    def __init__(self,name,start_last=False):
        self.name = name
        self.digital_outputs = {}
        self.digital_inputs = {}
        self.analogue_outputs = {}
        self.analogue_inputs = {}
        self.start_last=start_last
        inventory[self.name] = self
        
    def add_output(self,output,socket_number):
        if isinstance(output,DigitalOutput):
            self._add_to_outputdict(output,socket_number,self.digital_outputs,outtype='digital')
        elif isinstance(output,AnalogueOutput):
            self._add_to_outputdict(output,socket_number,self.analogue_outputs,outtype='analogue')
            
    def _add_to_outputdict(self,output,socket_number,outputdict,outtype):
            if socket_number >= self.n_outputs[outtype]:
                print 'WARNING:', self.type_string, self.name, 'only has', 
                print self.n_outputs[outtype], outtype, 'outputs.',
                print 'ignoring attempt to attach %s %s'%(output.type_string, output.name),
                print 'to', outtype,'output %s.\n'%str(socket_number)
                return
            if socket_number in outputdict.keys():
                existing_output = outputdict[socket_number]
                print 'WARNING:', self.type_string, self.name, 'already has', 
                print existing_output.type_string, existing_output.name,
                print 'attached to',outtype,'output %s.'%str(socket_number),
                print 'Overwriting with %s %s.\n'%(output.type_string, output.name)
            outputdict[socket_number] =  output
            
        
    def add_input(self,input,socket_number):
        pass
        
    def get_change_times(self):
        """Returns all the times at which on or more output values change."""
        # using a set ensures we get no duplicates:
        times = set()
        outputs = self.digital_outputs.values() + self.analogue_outputs.values()
        for output in outputs:
            for t in output.instructions.keys():
                times.add(t)
        times = list(times)
        return times     
       
    def create_timeseries(self):   
        times = self.get_change_times()  
        #let's find out what each output should be at each time:
        timeseries = {}
        for t in times:
            timeseries[t] = {}
        # digital outputs first:
        if self.digital_outputs:
            for t in times:
                digital_values = {}
                for socket_number, output in self.digital_outputs.items():
                    digital_values[socket_number] = output.get_value(t)
                timeseries[t]['digital'] = digital_values
        # Now for the analogue outputs:
        if self.analogue_outputs:
            for t in times:
                analogue_values = {}
                for socket_number, output in self.analogue_outputs.items():
                    analogue_values[socket_number] = output.get_value(t)
                timeseries[t]['digital'] = analogue_values
        return timeseries
    
    def program(self):
        """to be implemented by subclasses for specific devices"""
        pass
 
    def start(self):
        """to be implemented by subclasses for specific devices"""
        pass

    def stop(self):
        """to be implemented by subclasses for specific devices"""
        pass
 
class Output():
    type_string = 'generic output'
    def __init__(self,name,breakout_box,socket_number):
        self.name = name
        self.instructions = {}
        breakout_box.add_output(self,socket_number)
        
    def state_string(self,instructionvalue):
        return str(instructionvalue)
    
    def set_state(self,state,t):
        if t in self.instructions.keys():
            print 'WARNING: State of', self.type_string, self.name, 'at t=%ss'%str(t),
            print 'has already been set to %s.'%self.state_string(self.instructions[t]),
            print 'Overwriting to %s.\n'%self.state_string(state)
        self.instructions[t] = state
        
    def get_value(self,t):
        times = self.instructions.keys()
        # do we have an instruction at this exact time?
        if t in times:
            # if so, return the state for that instruction:
            return self.instructions[t]
        # if not, what was the preceding instruction?
        preceeding = max([time for time in times if time < t])
        return self.instructions[preceeding]
 
class AnalogueOutput(Output):
    type_string = 'analogue output'
    
        
class DigitalOutput(Output):
    type_string = 'digital output'
    def state_string(self,instructionvalue):
        return {1:'high',0:'low'}[instructionvalue]
  
    
class PulseBlaster(IODevice):
    import spinapi
    type_string = 'pulse blaster'
    n_outputs = {'analogue':2,'digital':12}
    n_inputs = {'analogue':None,'digital':None}
    
    def program(self):
        spinapi = self.spinapi
        if self.analogue_outputs:
            raise NotImplementedError('Sorry, only digital outputs for pulseblasters at the moment')
        timeseries = self.create_timeseries()
        times = timeseries.keys()
        print 'sorting'
        times.sort()
        print 'done sorting'
        spinapi.pb_init()
        spinapi.pb_core_clock(50)
        spinapi.pb_start_programming(spinapi.PULSE_PROGRAM)
        for j,t in enumerate(times):
            flags = ['0']*self.n_outputs['digital']
            for i in range(len(flags)):
                if i in timeseries[t]['digital'].keys():
                    flags[i] = str(timeseries[t]['digital'][i])
                    
            flagstring = ''.join(flags)
            if j > 0:
                preceding_time = times[j-1]
            else:
                preceding_time = -1 
            spinapi.pb_inst(flagstring,spinapi.CONTINUE,0,(t-preceding_time)*1e9)
        n_instructions = spinapi.pb_inst('000000000000',spinapi.STOP,0,1000)
        spinapi.pb_stop_programming()
        spinapi.pb_close()
        print 'programmed pulseblaster with',n_instructions,'instructions.'
        
    def start(self):
        spinapi = self.spinapi
        spinapi.pb_init()
        spinapi.pb_core_clock(50)
        spinapi.pb_start()
        
    def stop(self):
        spinapi = self.spinapi
        spinapi.pb_stop()
        spinapi.pb_close()

        
class Shutter(DigitalOutput):
    type_string = 'shutter'
    def state_string(self,instructionvalue):
        return {1:'open',0:'closed'}[instructionvalue]
    def open(self,t):
        self.set_state(1,t)
    def close(self,t):
        self.set_state(0,t)

def plot_timeseries(digital=True, analogue=True, IO_device_names=None):
    if IO_device_names is None:
        IO_device_names = inventory.keys()
    from pylab import plot, show,linspace, axis
    outputs = []
    times = []
    for device_name in IO_device_names:
        device = inventory[device_name]
        if digital:
            outputs.extend(device.digital_outputs.values())
        if analogue:
            outputs.extend(device.analogue_outputs.values())
        times.extend(device.get_change_times())
        
    times = linspace(min(times),max(times),len(times)*10)
    for output in outputs:
        plot(times,[output.get_value(t) for t in times])
    if outputs:
        axis([min(times),max(times),0,4])
        show()

def run():
    print 'running'
    last_device_to_start = []
    for device in inventory.values():
        device.program()
    for device in inventory.values():
        if device.start_last:
            last_device_to_start.append(device)
        else:
            device.start()
        if len(last_device_to_start) > 1:
            print 'WARNING: Multiple IO devices set to be started last:'
            for device in last_device_to_start:
                print '    %s %s'%(device.type_string,device.name)
            print ' %s %s'%(last_device_to-start[0].type_string,last_device_to-start[0].name),
            print 'will be the last to start running.'
            last_device_to_start[0].start()
            
def stop():
    for device in inventory.values():
        device.stop()
inventory = {}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
