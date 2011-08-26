from __future__ import division
import h5py
import sys
import itertools
import types

from pylab import *


def make_single_hdf_file(filename,argdict):
    print argdict
    f = h5py.File(filename,'w')
    group = f.create_group('params')
    for var in argdict:
        group.attrs[var] = argdict[var]
    f.close()
  
    
def make_hdf_files(basename, argdict):
    names = []  
    vals = []
    for key in argdict:
        if isinstance(argdict[key],types.GeneratorType):
           result = list(argdict[key])
        elif isinstance(argdict[key],tuple):
           start,stop,n = argdict[key]
           result = linspace(start,stop,n)
        elif isinstance(argdict[key], ndarray) or  isinstance(argdict[key], list):
            result = argdict[key]
        else:
            result = [argdict[key]]
        names.append(key)
        vals.append(result)
  
    nruns = 1
    for lst in vals:
        nruns *= len(lst)
    ndigits = int(ceil(log10(nruns)))
    print ndigits
    print nruns
    for i, values in enumerate(itertools.product(*vals)):
        rundict = {} 
        for name,val in zip(names,values):
            rundict[name] = val
        make_single_hdf_file(('%s%0'+str(ndigits)+'d.h5')%(basename,i),rundict)
            
if __name__ == '__main__':
    try:
        assert len(sys.argv) > 1
        basename = sys.argv[-1]
        assert basename.lower().endswith('h5')
    except:
        sys.stderr.write('ERROR: No output hdf5 filename provided. Stopping.\n')
        sys.exit(1)
      
    argdict = {}     
    for arg in sys.argv[1:-1]:
        key,val = arg.split('=')
        argdict[key] = eval(val)
    make_hdf_files(basename.split('.h5')[0], argdict)
      
