from labscript import *

# Initialisation, runs at import. Can be suppressed by setting
# labscript_auto_init = False in the locals of the importing scope
# before importing labscript. If you do this, you'll need to call
# labscript_init() yourself:

import inspect
importing_frame = inspect.currentframe()
importing_locals = importing_frame.f_back.f_locals
if not 'labscript_auto_init' in importing_locals or importing_locals['labscript_auto_init']:
    print 'auto init!'
    if len(sys.argv) > 1:
        labscript_init(sys.argv[1],labscript_file=sys.argv[0])
    elif sys.argv[0]:
        labscript_init(sys.argv[0].replace('.py','.h5'), labscript_file=sys.argv[0], new=True)
