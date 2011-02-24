import analoguetest
import time

class NIBoard:
    description = 'generic NI board'
    clock_limit = 1e6 #TODO get a real number to put here
    def __init__(self, name, clockedby):
        self.name = name
        if clockedby.clock_limit > self.clock_limit:
            clockedby.clock_limit = self.clock_limit
        self.outputs = clockedby.outputs
       
class PulseBlaster(analoguetest.IODevice):
    description = 'PulseBlaster'
    clock_limit = 25.0e6 # TODO get a real number to put here
    clock_connection = 11
    direct_outputs = []
    def generate_code(self):
        self.make_instruction_table()
        for output in self.outputs:
            if output.connected_to_device is self:
                self.direct_outputs.append(output)
               
               
inventory = []
if __name__ == '__main__':
    start_time = time.time()

    pulseblaster1 = PulseBlaster('PulseBlaster_1')
    NI_board1 = NIBoard('NI_board_1', pulseblaster1)
    output1 = analoguetest.Output('output 1',NI_board1,0)
    output2 = analoguetest.Output('output 2',NI_board1,1)
    output3 = analoguetest.Output('output 3',NI_board1,2)

    output1.add_instruction(0,2)
    output1.add_instruction(1, {'function': analoguetest.ramp(1,2,2,3), 'end time' : 3, 'clock rate':  5})

    output2.add_instruction(0,3)
    output2.add_instruction(2, {'function': analoguetest.ramp(2,3,3,4), 'end time' : 5, 'clock rate': 10})
    output2.add_instruction(5.9,5)
    output2.add_instruction(7,4)
    output2.add_instruction(8,5)
    output3.add_instruction(0, {'function': analoguetest.sine(0,1), 'end time' : 10, 'clock rate':    100000})

    pulseblaster1.generate_code()
    print time.time() - start_time
    analoguetest.plot_outputs()
