import h5py
from pylab import *
import sys

f = h5py.File('dummy.h5','w')
group = f.create_group('params')
group.attrs['z'] = 'test'
group.attrs['A'] = 7
group.attrs['B'] = 1e6
group.attrs['C'] = 'hello!'
group.attrs['D'] = True
#group.attrs['lambda'] = 1
f.close()
