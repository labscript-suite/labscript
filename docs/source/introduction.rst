Introduction
=====================================
The labscript API is used to define the logic of an experiment that you wish to run. It is recommended that you read our paper before this documentation, so you are familiar with terms like pseudoclock. It would also be a good idea to familiarise yourself with the Python programming language and object oriented (OO) programming if you are not already. 

To give you an idea of what a sample experiment looks like, the simplest experiment script (that does something) using the labscript API is below::

    from labscript import *
    from labscript_devices.PulseBlaster import PulseBlaster
    
    # Connection Table
    PulseBlaster(name='pulseblaster_0', board_number=0)
    DigitalOut(name='my_digital_out', parent_device=pulseblaster_0.direct_outputs, connection='flag 2')

    #Experiment Logic
    start()
    my_digital_out.go_low(t=0)  # start low at the start
    my_digital_out.go_high(t=1) # go high at 1s
    stop(2)                     # stop at 2s
    
The script consists of two parts, the connection table and the experiment logic which will be discussed in the following sections.


    
