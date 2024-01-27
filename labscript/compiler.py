#####################################################################
#                                                                   #
# /compiler.py                                                      #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

"""The labscript compiler interface. This is only relevant to developers and those
interested in the labscript interface to runmanager."""

import builtins

from labscript_utils.labconfig import LabConfig

_builtins_dict = builtins.__dict__

# Extract settings from labconfig
_SAVE_HG_INFO = LabConfig().getboolean("labscript", "save_hg_info", fallback=False)
_SAVE_GIT_INFO = LabConfig().getboolean("labscript", "save_git_info", fallback=False)


class Compiler(object):
    """Compiler singleton that saves relevant parameters during compilation of each shot"""

    _instance = None

    _existing_builtins_dict = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            # Save initial builtins state so that we (and possibly someone else) can
            # call ``reset()`` at any time.
            cls._instance.save_builtins_state()
            # Initialise state
            cls._instance.reset()
        return cls._instance

    def save_builtins_state(self):
        self._existing_builtins_dict = _builtins_dict.copy()

    def reset(self):
        """restores builtins and the labscript module to its state before
        labscript_init() was called"""
        # Reset builtins
        for name in _builtins_dict.copy(): 
            if name not in self._existing_builtins_dict:
                del _builtins_dict[name]
            else:
                _builtins_dict[name] = self._existing_builtins_dict[name]

        # Reset other variables

        # The labscript file being compiled:
        self.labscript_file = None
        # All defined devices:
        self.inventory = []
        # The filepath of the h5 file containing globals and which will
        # contain compilation output:
        self.hdf5_filename = None

        self.start_called = False
        self.wait_table = {}
        self.wait_monitor = None
        self.master_pseudoclock = None
        self.all_pseudoclocks = None
        self.trigger_duration = 0
        self.wait_delay = 0
        self.time_markers = {}
        self._PrimaryBLACS = None
        self.save_hg_info = _SAVE_HG_INFO
        self.save_git_info = _SAVE_GIT_INFO
        self.shot_properties = {}

        # This used to be in a separate config object, but it's been moved here so it
        # gets reset
        self.suppress_mild_warnings = True
        self.suppress_all_warnings = False
        self.compression = 'gzip'  # set to 'gzip' for compression 


compiler = Compiler()
"""The compiler instance"""
