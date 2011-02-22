from pylab import *

def ramp(t,initial,duration,final):
    m = (final - initial)/float(duration) # be sure to prevent integer division!
    c = initial - m*t
    return lambda x: m*x + c

def sine(t,frequency):
    return lambda x: sin(x) + 0.5

class IODevice:
    def __init__(self,name):
        self.outputs = []
        self.name = name
        
    def collect_change_times(self):
        # use a set to avoid duplicates:
        self.change_times = set()
        for output in self.outputs:
            output.make_times()
            self.change_times = self.change_times.union(output.times)
        self.change_times = list(self.change_times)
        self.change_times.sort()
        
    def make_timeseries(self):
        for output in self.outputs:
            output.make_timeseries(self.change_times)
        
    def expand_change_times(self):
        self.all_times = []
        for i, time in enumerate(self.change_times):
            # what's the fastest clock rate?
            maxrate = 0
            for output in self.outputs:
                # check if output is sweeping and has highest clock rate so far:
                if isinstance(output.timeseries[i],dict) and output.timeseries[i]['clock rate'] > maxrate:
                    maxrate = output.timeseries[i]['clock rate']
            if maxrate:
                # If there was sweeping at this timestep, store an array of times at the max clock rate:
                n_ticks = int((self.change_times[i+1] - time)*maxrate)
                duration = n_ticks/float(maxrate) # avoiding integer division
                self.all_times.append(linspace(time,time + duration,n_ticks,endpoint=False))
            else:
                self.all_times.append(time)
                
    def expand_timeseries(self):
        for output in self.outputs:
            output.make_outputarray(self.all_times)
            
    def flatten(self,inarray,total_points):
        flat = zeros(total_points)
        i = 0
        for val in inarray:
            if iterable(val):
                flat[i:i+len(val)] = val[:]
                i += len(val)
            else:
                flat[i] = val
                i += 1
        return flat
                
    def make_raw_output(self):
        total_points = sum([len(times) if iterable(times) else 1 for times in self.all_times])
        self.flat_times = self.flatten(self.all_times,total_points)
        for output in self.outputs:
            output.raw_output = self.flatten(output.outputarray,total_points)
        
    def make_instruction_table(self):
        self.collect_change_times()
        self.make_timeseries()
        self.expand_change_times()
        self.expand_timeseries()
        self.make_raw_output()
        
class Output:
    def __init__(self,name,IO_device,outputnumber):
        self.name = name
        self.instructions = {}
        IO_device.outputs.append(self)
        
    def add_instruction(self,time,instruction):
        self.instructions[time] = instruction
 
    def make_times(self):
        # Check that ramps have instructions following them.
        # If they don't, insert an instruction telling them to hold their final value.
        for instruction in self.instructions.values():
            if isinstance(instruction, dict) and instruction['end time'] not in self.instructions.keys():
                self.add_instruction(instruction['end time'], instruction['function'](instruction['end time']))
        # Check if there are no instructions. Generate a warning and insert an
        # instruction telling the output to remain at zero.
        if not self.instructions:
            print 'WARNING:', self.name, 'has no instructions. It will be set to output zero for all time.'
            self.add_instruction(0,0)
        self.times = self.instructions.keys()
        self.times.sort()
        
    def make_timeseries(self,change_times):
        self.timeseries = []
        i = 0
        for change_time in change_times:
            try:
                if i < len(self.times):
                    while change_time >= self.times[i]:
                        i += 1
            except IndexError:
                pass    
            instruction = self.instructions[self.times[i-1]]
            if isinstance(instruction, dict) and instruction['end time'] <= change_time:
                print 'detected that a ramp has ended!' 
                instruction = instruction['function'](instruction['end time'])
            self.timeseries.append(instruction)     
    
    def make_outputarray(self,all_times):
        self.outputarray = []
        for i, time in enumerate(all_times):
            if iterable(time):
                if isinstance(self.timeseries[i],dict):
                    outarray = self.timeseries[i]['function'](time)
                else:
                    outarray = self.timeseries[i]*ones(len(time))
                self.outputarray.append(outarray)
            else:
                self.outputarray.append(self.timeseries[i])
    
    
def main():
    import time
    start_time = time.time()
    
    device1 = IODevice('device 1')

    output1 = Output('output 1',device1,1)
    output2 = Output('output 2',device1,2)
    output3 = Output('output 3',device1,3)

    output1.add_instruction(0,2)
    output1.add_instruction(1, {'function': ramp(1,2,2,3), 'end time' : 3, 'clock rate': 5})
    #output1.add_instruction(3,3)

    output2.add_instruction(0,3)
    output2.add_instruction(2, {'function': ramp(2,3,3,4), 'end time' : 5, 'clock rate': 10})
    output2.add_instruction(5,4)
    output2.add_instruction(6,5)
    output2.add_instruction(7,4)
    output2.add_instruction(8,5)
    

    device1.make_instruction_table()
    print time.time() - start_time
    pointstrings = ['ro-','bo-','go-']
    for i, output in enumerate(device1.outputs):
        plot(device1.flat_times,output.raw_output,pointstrings[i],label=output.name)
    grid(True)
    xlabel('time (seconds)')
    ylabel('analogue output values')
    title('Putting analogue outputs on a common clock')
    legend(loc='lower right')
    axis([0,10,-1,5.5])
    show()
#    



main()
