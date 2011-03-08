import time
start_time = time.time()    
import sys
import h5py
import socket  
host='130.194.171.157'
port=6000

if not len(sys.argv) > 1:
    sys.stderr.write('ERROR: No hdf5 file provided as a command line argument. Stopping.\n')
    sys.exit(1)
    
try:
    novatech = socket.create_connection((host,port), timeout=1)
except:
    sys.stderr.write('ERROR: could not connect to Novatech DDS on %s:%s. Stopping.\n'%(host,port))
    sys.exit(1)


with h5py.File(sys.argv[-1],'r') as hdf5_file:
    instructions = hdf5_file['/Novatech DDS/TABLE_DATA']
    novatech.send('m 0\n')
    for i, (freq0,freq1,phase0,phase1,amp0,amp1) in enumerate(instructions[:len(instructions)-1]):
        novatech.send('t0 %04x %08x,%04x,%04x,ff\n'%(i,freq0,phase0,amp0))
        novatech.send('t1 %04x %08x,%04x,%04x,ff\n'%(i,freq1,phase1,amp1))
    i, (freq0,freq1,phase0,phase1,amp0,amp1) = len(instructions) - 1, instructions[-1]
    novatech.send('t0 %04x %08x,%04x,%04x,00\n'%(i,freq0,phase0,amp0))
    novatech.send('t1 %04x %08x,%04x,%04x,00\n'%(i,freq1,phase1,amp1))
    novatech.send('m t\n')
    novatech.close()
print 'programming Novatech DDS:',round(time.time() - start_time,2),'sec'
