import sys
import h5py
hdf5_filename = sys.argv[1]
outfile = sys.argv[2]

hdf5_file = hdf5_file = h5py.File(hdf5_filename)
connection_table = hdf5_file['/connection table']
field = 'name'
fieldlengths = []
for field in ['name','class','parent','connected to']:
    fieldlengths.append(max([len(item) for item in connection_table[field]]))
lines = ['from labscript import *','',]
for name,class_,parent,connected_to in connection_table:
   lines.append('%s('%class_.ljust(fieldlengths[1] + 1) + \
          ('\'%s\', '%name).rjust(fieldlengths[0] + 5) + \
          ('%s, '%parent).rjust(fieldlengths[2] + 2) + \
          ('\'%s\' )'%connected_to if connected_to != 'None' else 'None )').rjust(fieldlengths[3]+4))
          
lines.extend(['','stop(t=0)'])
script =  '\n'.join(lines)

open(outfile,'w').write(script)
