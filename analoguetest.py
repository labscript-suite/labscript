from pylab import *

def ramp(starttime,startvalue,duration,endvalue):
    m = (endvalue - startvalue)/float(duration) # be sure to prevent integer division!
    c = startvalue - m*starttime
    return lambda t: m*t + c

class IODevice:
    def __init__(self):
        self.outputs = []
        
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
                if isinstance(output.timeseries[i],list) and output.timeseries[i][1] > maxrate:
                    maxrate = output.timeseries[i][1]
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
        for output in self.outputs:
            plot(self.flat_times,output.raw_output,'ro')
        show()
        
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
            self.timeseries.append(self.instructions[self.times[i-1]])     
    
    def make_outputarray(self,all_times):
        self.outputarray = []
        for i, time in enumerate(all_times):
            if iterable(time):
                if iterable(self.timeseries[i]):
                    outarray = self.timeseries[i][0](time)
                else:
                    outarray = self.timeseries[i]*ones(len(time))
                self.outputarray.append(outarray)
            else:
                self.outputarray.append(self.timeseries[i])
    
    def make_raw_output(self,total_points):
        raw_output = zeros(total_points)
        i = 0
        for outvalue in self.outputarray:
            if iterable(outvalue):
                raw_output[i:i+len(outvalue)] = outvalue[:]
                i += len(outvalue)
            else:
                raw_output[i] = outvalue
                i += 1
        print raw_output
    
def main():
    device1 = IODevice()

    output1 = Output('output1',device1,1)
    output2 = Output('output2',device1,2)

    output1.add_instruction(0,2)
    output1.add_instruction(1,[ramp(1,2,2,3),10])
    output1.add_instruction(3,3)

    output2.add_instruction(0,3)
    output2.add_instruction(2,[ramp(2,3,2,4),20])
    output2.add_instruction(5,4)
    output2.add_instruction(3,3.5)
    output2.add_instruction(4,4)
    output2.add_instruction(5,3.5)
    output2.add_instruction(6,4)
    output2.add_instruction(7,3.5)

    device1.make_instruction_table()




main()





#class Test:
#    def make_timeseries(self,change_times):
#        self.timeseries = []
#        i = 0
#        for change_time in change_times:
#            if i < len(self.times) - 1:
#                while change_time >= self.times[i]:
#                    i += 1
#            self.timeseries.append(self.instructions[self.times[i-1]])     
#        print self.timeseries
# 
#test = Test()       
#change_times = [1,2,3,4,5,6,7,8,9]
#test.times = [1,3,5,7]
#test.instructions = {1:'one',3:'three',5:'five',7:'seven'}
#test.timeseries = []
#test.make_timeseries(change_times)

#i = 0
#for change_time in change_times:
#    print i
#    print times[i]
#    if i < len(times)-1:
#        while change_time >= times[i]:
#            i += 1
#    timeseries.append(instructions[times[i-1]])    
#    print timeseries
        
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
        
#from pylab import *
#myramp = ramp(starttime=1,startvalue=2,duration=3,endvalue=4)
#t = linspace(1,4,100)
#plot(t,myramp(t))
