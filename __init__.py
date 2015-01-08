#####################################################################
#                                                                   #
# /__init__.py                                                      #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

from labscript import *

try:
    from labscript_utils import check_version
except ImportError:
    raise ImportError('Require labscript_utils > 2.1.0')

check_version('labscript_utils', '2.2', '3')


# Initialisation, runs at import. Can be suppressed by setting
# labscript_auto_init = False in the locals of the importing scope
# before importing labscript. If you do this, you'll need to call
# labscript_init() yourself:

#import inspect
#importing_frame = inspect.currentframe()
#importing_locals = importing_frame.f_back.f_locals
#if not 'labscript_auto_init' in importing_locals or importing_locals['labscript_auto_init']:
#    overwrite = False
#    if '-f' in sys.argv:
#        overwrite = True
#        sys.argv.remove('-f')
#    if len(sys.argv) > 1:
#        labscript_init(sys.argv[1],labscript_file=sys.argv[0])
#    elif sys.argv[0]:
#        labscript_init(sys.argv[0].replace('.py','.h5'), labscript_file=sys.argv[0], new=True, overwrite=overwrite)

__version__ = '2.1.0'
