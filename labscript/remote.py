#####################################################################
#                                                                   #
# /remote.py                                                        #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

"""Classes for configuring remote/secondary BLACS and/or device workers"""

from .compiler import compiler
from .labscript import Device, set_passed_properties


class _PrimaryBLACS(Device):
    pass


class _RemoteConnection(Device):
    @set_passed_properties(
        property_names = {}
    )
    def __init__(self, name, parent=None, connection=None):
        if parent is None:
            # define a hidden parent of top level remote connections so that
            # "connection" is stored correctly in the connection table
            if compiler._PrimaryBLACS is None:
                compiler._PrimaryBLACS = _PrimaryBLACS('__PrimaryBLACS', None, None)
            parent = compiler._PrimaryBLACS
        Device.__init__(self, name, parent, connection)


class RemoteBLACS(_RemoteConnection):
    def __init__(self, name, host, port=7341, parent=None):
        _RemoteConnection.__init__(self, name, parent, "%s:%s"%(host, port))


class SecondaryControlSystem(_RemoteConnection):
    def __init__(self, name, host, port, parent=None):
        _RemoteConnection.__init__(self, name, parent, "%s:%s"%(host, port))
