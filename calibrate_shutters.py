# This file just creates a dummy calibration file for the moment. This
# will get more sophisticated as time goes on. Ultimately there will be
# a calibration experiment script, and an analysis script. The experiment
# script will open and close all the shutters and take power measurements
# from photodiodes in the laser paths. The analysis script will take the
# resulting hdf5 file as input and fit the laser power measurements to
# step functions, comparing the fitted step time to the time the shutter
# was told to open. Then it will save the results to a hdf5 file. If
# there was an existing calibration file, it will replace only the
# calibrations that were being redone, keeping calibrations for other
# devices. The old calibration file will be copied to a folder where it
# can be kept permanently. That's the idea anyway.

import h5py
f = h5py.File('calibrations.h5','w')
g = f.create_group('shutters')
shutter1 = g.create_group('shutter 1')
shutter2 = g.create_group('shutter 2')
shutter1.attrs['open_delay'] = 5e-3
shutter1.attrs['close_delay'] = 5e-3
shutter2.attrs['open_delay'] = 7.5e-3
shutter2.attrs['close_delay'] = 7.5e-3
f.close()
