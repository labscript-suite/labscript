usage = """USAGE:  python novatech.py [-serial=com1 | -tcp=192.168.1.10:10001] [-program] [-verify] -boardnumber=1 infile.h5"""

import time   
import sys
import h5py
import serial
import socket
import os

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
           
def send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell):
    novatech.write('t0 %04x %08x,%04x,%04x,%s\r\n'%(i,freq0,phase0,amp0,dwell))
    print i, repr(novatech.readline()),
    novatech.write('t1 %04x %08x,%04x,%04x,%s\r\n'%(i,freq1,phase1,amp1,dwell))
    print i, repr(novatech.readline())

def verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell) :
    success = True
    novatech.write('d0 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell):
        print  'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell)
        success = False
    else:
        print  'SUCCESS',
    novatech.write('d1 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell):
        print  'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell)
        success = False
    else:
        print  'SUCCESS'
    return success

if not ('-program' in sys.argv or '-verify' in sys.argv):
    sys.stderr.write('novatech.py has not been instructed to program or verify. Stopping.\n')
    print usage  
    sys.exit(1)

for i, arg in enumerate(sys.argv):
    if arg.startswith('-boardnumber'):
        boardnumber = arg.split('=')[1]
        break
    if i == len(sys.argv) - 1:
        sys.stderr.write('ERROR: No board number provided. Stopping.\n')   
        print usage
        sys.exit(1)
            
if not len(sys.argv) > 4:
    sys.stderr.write('ERROR: Require hdf5 file. Stopping.\n')
    print usage
    sys.exit(1)

if not os.path.exists(sys.argv[-1]):
    sys.stderr.write('hdf5 file %s does not exist. Stopping.\n'%sys.argv[-1])
    print usage
    sys.exit(1)
    
for i, arg in enumerate(sys.argv):  
    if arg.startswith('-serial'):
        method = 'serial'
        try:
            com = int(arg.lower().replace('-serial=com','')) - 1
        except:
            sys.stderr.write('ERROR: problem parsing command line arg for serial com port. Stopping.\n')
            print usage
            raise
        novatech = connect_serial()
        break
    elif arg.startswith('-tcp'):
        method = 'serial'
        try:
            host,port = arg.lower().replace('-tcp=','').split(':')
            port = int(port)
        except:
            sys.stderr.write('ERROR: problem parsing command line arg for tcp/ip address. Stopping.\n')
            print usage
            raise
        novatech = connect_tcp()
        break
    if i == len(sys.argv) - 1:
        sys.stderr.write('ERROR: No connection method (serial or tcp/ip) provided. Stopping.\n')   
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
            novatech.close()
            sys.exit(1)
        continue

success = True               
with h5py.File(sys.argv[-1],'r') as hdf5_file:
    instructions = hdf5_file['/NT-DDS9M_%s/TABLE_DATA'%boardnumber]
    if '-program' in sys.argv:
        for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions):
            send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00' if i == len(instructions)-1 else 'ff')
        print 'programming Novatech DDS:', round(time.time() - start_time,2),'sec'
    if '-verify' in sys.argv:
        start_time = time.time()
        for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions):
            if not verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00' if i == len(instructions)-1 else 'ff'):
                success = False
        print 'Verifying instructions of Novatech DDS:', round(time.time() - start_time,2),'sec'
if success == False:
    sys.stderr.write('ERROR: One or more instructions were not verified to have been programmed correctly. Stopping.\n')
    novatech.close()
    sys.exit(1)
    
novatech.close()



