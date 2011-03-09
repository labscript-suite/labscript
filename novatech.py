import time
start_time = time.time()    
import sys
import h5py
import socket  
host='130.194.171.157'
port=10001
wait_time = 0.03
if not len(sys.argv) > 1:
    sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
    sys.exit(1)
class dummy:
    def send(self,string):
        print string,
try:
    novatech = socket.create_connection((host,port), timeout=1)
except:
    sys.stderr.write('ERROR: could not connect to Novatech DDS on %s:%s. Stopping.\n'%(host,port))
    sys.exit(1)

responding = False
while not responding:
    novatech.send('e d\n')
    try:
        response = novatech.recv(512).strip()
    except:
        continue
    print 'e d', response
    if response == 'OK':
        responding = True
novatech.send('e d\n')
print 'i e', novatech.recv(8).strip()
with h5py.File(sys.argv[-1],'r') as hdf5_file:
    instructions = hdf5_file['/Novatech DDS/TABLE_DATA']
    novatech.send('m 0\n')
    print 'm 0', novatech.recv(64).strip()
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        novatech.send('t0 %04x %08x,%04x,%04x,ff\n'%(i,freq0,phase0,amp0))
        time.sleep(wait_time)
        novatech.send('t1 %04x %08x,%04x,%04x,ff\n'%(i,freq1,phase1,amp1))
        print i#, novatech.recv(8).strip()
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    novatech.send('t0 %04x %08x,%04x,%04x,00\n'%(i,freq0,phase0,amp0))
    time.sleep(wait_time)
    novatech.send('t1 %04x %08x,%04x,%04x,00\n'%(i,freq1,phase1,amp1))
    time.sleep(wait_time)
    novatech.send('m t\n')
    time.sleep(wait_time)
    print 'm t', novatech.recv(4096).replace('\r\n','')
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        novatech.send('d0 %04x\n'%i)
        time.sleep(wait_time)
        result = novatech.recv(512).strip().replace(' OK','')
        if not result.lower() == '%08x,%04x,%04x,ff'%(freq0,phase0,amp0):
            print i, 'FAIL:', result.lower(), '|| %08x,%04x,%04x,ff'%(freq0,phase0,amp0)
        else:
            print i, 'SUCCESS'
        novatech.send('d1 %04x\n'%i)
        time.sleep(wait_time)
        result = novatech.recv(512).strip().replace(' OK','')
        if not result.lower() == '%08x,%04x,%04x,ff'%(freq1,phase1,amp1):
            print i, 'FAIL:', result.lower(), '|| %08x,%04x,%04x,ff'%(freq1,phase1,amp1)
        else:
            print i, 'SUCCESS'
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    novatech.send('d0 %04x\n'%i)
    time.sleep(wait_time)
    result = novatech.recv(512).strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,00'%(freq0,phase0,amp0):
        print i, 'FAIL:', result.lower(), '|| %08x,%04x,%04x,00'%(freq0,phase0,amp0)
    else:
        print i, 'SUCCESS'
    novatech.send('d1 %04x\n'%i)
    time.sleep(wait_time)
    result = novatech.recv(512).strip().replace(' OK','')
    if not result.lower() == '%08x,%04x,%04x,00'%(freq1,phase1,amp1):
        print i, 'FAIL:', result.lower(), '|| %08x,%04x,%04x,00'%(freq1,phase1,amp1)
    else:
        print i, 'SUCCESS'
    novatech.close()
print 'programming Novatech DDS:',round(time.time() - start_time,2),'sec'













