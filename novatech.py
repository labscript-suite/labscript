usage = """USAGE:  python novatech.py [-serial=com1 | -tcp=192.168.1.10:10001] [-program] [-verify] infile.h5"""

import time   
import sys
import h5py
import serial
import socket

class SocketWithSerialMethods(object):
    def __init__(self,addr, timeout=1):
        self.s = socket.create_connection(addr, timeout=timeout)
    def write(self,text):
        self.s.send(text)
    def readline(self, n=512):
        try: 
            return self.s.recv(n)
        except KeyboardInterrupt:
            raise
        except:
            return ''
    def close(self):
        self.s.close()

def connect_serial():
    try:
        novatech = serial.Serial(com, baudrate = 115200, timeout=1)
    except:
        sys.stderr.write('couldn\'t connect to Novatech DDS on COM %s. Stopping.\n'%str(com+1))
        sys.exit(0)
    return novatech
    
def connect_tcp():
    try:
        novatech = SocketWithSerialMethods((host,port),timeout=2)
    except:
        raise
        sys.stderr.write('couldn\'t connect to Novatech DDS on %s:%s. Stopping.\n'%(host,str(port)))
        sys.exit(0)
    return novatech
           
if not len(sys.argv) > 2:
    sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
    print usage
    sys.exit(1)


if sys.argv[1].startswith('-serial'):
    method = 'serial'
    try:
        com = int(sys.argv[1].lower.replace('-serial=com','')) - 1
    except:
        print usage
        raise
    novatech = connect_serial()
elif sys.argv[1].startswith('-tcp'):
    method = 'serial'
    try:
        host,port = sys.argv[1].lower().replace('-tcp=','').split(':')
        port = int(port)
        print host,port
    except:
        print usage
        raise
    novatech = connect_tcp()
else:
    print usage
    sys.exit(1)
    
start_time = time.time()  
responding = False
i = 0
while not responding:
    novatech.write('e d\n')
    response = novatech.readline()
    print 'e d', repr(response)
    if response == 'e d':
        response == novatech.readline()
    if response == 'OK\r\n':
        responding = True
        print 'NovaTech DDS is responding properly.'
    else:
        i += 1
        if i == 10:
            sys.stderr.write('NovaTech DDS not responding to commands. Stopping.\n')
            sys.exit(1)
        continue
            
def send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell):
    novatech.write('t0 %04x %08x,%04x,%04x,%s\r\n'%(i,freq0,phase0,amp0,dwell))
    print i, repr(novatech.readline()),
    novatech.write('t1 %04x %08x,%04x,%04x,%s\r\n'%(i,freq1,phase1,amp1,dwell))
    print i, repr(novatech.readline())

def verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell) :
    novatech.write('d0 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell):
        print  'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell)
    else:
        print  'SUCCESS',
    novatech.write('d1 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell):
        print  'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell)
    else:
        print  'SUCCESS'
if not ('-program' in sys.argv or '-verify' in sys.argv):
    print 'novatech.py has no instructions to program or verify. Stopping.'
    print usage  
    novatech.close()
    sys.exit(1)
    
with h5py.File(sys.argv[-1],'r') as hdf5_file:
    instructions = hdf5_file['/Novatech DDS/TABLE_DATA']
    if '-program' in sys.argv:
        for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions):
            send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00' if i == len(instructions)-1 else 'ff')
        print 'programming Novatech DDS:',round(time.time() - start_time,2),'sec'
    if '-verify' in sys.argv:
        for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions):
            verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00' if i == len(instructions)-1 else 'ff')

novatech.close()
