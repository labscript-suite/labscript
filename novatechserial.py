import time   
import sys
import h5py
import serial
import socket
host='130.194.171.157'
port=10001
com = 0
timeout = 2
wait_time = 0.05
class SocketWithSerialMethods(object):
    def __init__(self):
        self.s = socket.create_connection((host,port), timeout=timeout)
    def write(self,text):
        self.s.send(text)
    def readline(self, num=512):
        return self.s.recv(num)
    def close(self):
        self.s.close()
        
if '-serial' in sys.argv:
    try:
        novatech = serial.Serial(com, baudrate = 115200, timeout=timeout)
    except:
        sys.stderr.write('couldn\'t connect to Novatech DDS on COM %s. Stopping.\n'%str(com+1))
        sys.exit(0)
else:
    i = 0
    connected = False
    while not connected:
        try:
            novatech = SocketWithSerialMethods()
            connected = True
        except:
            i += 1
            if i == 5:
                break
            continue
    if i == 5:
        sys.stderr.write('Couldn\'t connect to Novatech DDS on %s:%s. Stopping.\n'%(host,str(port)))
        sys.exit(0)
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

start_time = time.time()            
novatech.write('e d\n')
response = novatech.readline().strip()
print 'e d', response
if not response == 'OK':
    sys.stderr.write('NovaTech DDS not responding to commands. Stopping.\n')
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
    try:
        print 'results:', novatech.readline(4096*2)
    except:
        pass
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'ff')  
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    verify_an_instruction(i, freq0,freq1,phase0,phase1,amp0,amp1,'00')
print 'programming Novatech DDS:',round(time.time() - start_time,2),'sec'
novatech.close()
