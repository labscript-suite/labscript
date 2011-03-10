import time
start_time = time.time()    
import sys
import h5py
import serial
import socket
host='130.194.171.157'
port=10001
com = 3

class SocketWithSerialMethods(object):
    def __init__(self,(host,port),timeout):
        self.sock = socket.create_connection((host,port), timeout=1)
    def write(self,text):
        self.sock.send(text)
    def readline(self):
        return self.sock.recv(512)
        
if '-serial' in sys.argv:
    novatech = serial.Serial(com, baudrate = 115200, timeout=1)
else:
    novatech = SocketWithSerialMethods((host,port), timeout=2)
    responding = False
    while not responding:
        novatech.write('e d\n')
        try:
            response = novatech.readline().strip()
        except:
            continue
        print 'e d', response
        if response == 'OK':
            responding = True

            
novatech.write('e d\n')
response = novatech.readline().strip()
print 'e d', response
if not response == 'OK':
    sys.syderr.write('NovaTech DDS not responding to commands. Stopping.\n')
    sys.exit(1)

def send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell):
    novatech.write('t0 %04x %08x,%04x,%04x,%s\n'%(i,freq0,phase0,amp0,dwell))
    print i, novatech.readline().strip()
    novatech.write('t1 %04x %08x,%04x,%04x,%s\n'%(i,freq1,phase1,amp1,dwell))
    print i, novatech.readline().strip()

def verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,dwell) :
    novatech.write('d0 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell):
        print i, 'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq0,phase0,amp0,dwell)
    else:
        print i, 'SUCCESS'
    novatech.write('d1 %04x\n'%i)
    result = novatech.readline().strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell):
        print i, 'FAIL:', result.lower(), '!= %08x,%04x,%04x,%s'%(freq1,phase1,amp1,dwell)
    else:
        print i, 'SUCCESS'
        
with h5py.File(sys.argv[-1],'r') as hdf5_file:
    instructions = hdf5_file['/Novatech DDS/TABLE_DATA']
    novatech.write('m 0\n')
    print 'm 0', novatech.readline().strip()
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'ff')
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    send_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00')
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'ff')  
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00')
print 'programming Novatech DDS:',round(time.time() - start_time,2),'sec'
