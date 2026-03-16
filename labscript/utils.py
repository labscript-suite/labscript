#####################################################################
#                                                                   #
# /utils.py                                                         #
#                                                                   #
# Copyright 2013, Monash University                                 #
#                                                                   #
# This file is part of the program labscript, in the labscript      #
# suite (see http://labscriptsuite.org), and is licensed under the  #
# Simplified BSD License. See the license.txt file in the root of   #
# the project for the full license.                                 #
#                                                                   #
#####################################################################

"""Utility functions"""

import contextlib
from inspect import getcallargs
from functools import wraps

import numpy as np

from .compiler import compiler

cached_RemoteConnection = None
cached_ClockLine = None
cached_PseudoclockDevice = None


def print_time(t, description):
    """Print time with a descriptive string.

    Useful debug tool to print time at a specific point
    in the shot, during shot compilation. Helpful when
    the time is calculated.

    Args:
        t (float): Time to print
        description (str): Descriptive label to print with it
    """
    print(f"t={t:.9f}s: {description}")


def is_remote_connection(connection):
    """Returns whether the connection is an instance of ``_RemoteConnection``
    
    This function defers and caches the import of ``_RemoteConnection``. This both
    breaks the circular import between ``Device`` and ``_RemoteConnection``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    global cached_RemoteConnection
    if cached_RemoteConnection is None:
        from .remote import _RemoteConnection
        cached_RemoteConnection = _RemoteConnection
    return isinstance(connection, cached_RemoteConnection)


def is_clock_line(device):
    """Returns whether the connection is an instance of ``ClockLine``
    
    This function defers and caches the import of ``ClockLine``. This both
    breaks the circular import between ``Device`` and ``ClockLine``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    global cached_ClockLine
    if cached_ClockLine is None:
        from .core import ClockLine
        cached_ClockLine = ClockLine
    return isinstance(device, cached_ClockLine)


def is_pseudoclock_device(device):
    """Returns whether the connection is an instance of ``PseudoclockDevice``
    
    This function defers and caches the import of ``PseudoclockDevice``. This both
    breaks the circular import between ``Device`` and ``PseudoclockDevice``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    global cached_PseudoclockDevice
    if cached_PseudoclockDevice is None:
        from .core import PseudoclockDevice
        cached_PseudoclockDevice = PseudoclockDevice
    return isinstance(device, cached_PseudoclockDevice)


def set_passed_properties(property_names=None):
    """
    Decorator for device __init__ methods that saves the listed arguments/keyword
    arguments as properties. 

    Argument values as passed to __init__ will be saved, with
    the exception that if an instance attribute exists after __init__ has run that has
    the same name as an argument, the instance attribute will be saved instead of the
    argument value. This allows code within __init__ to process default arguments
    before they are saved.

    Internally, all properties are accessed by calling :obj:`self.get_property() <Device.get_property>`.
    
    Args:
        property_names (dict): A dictionary {key:val}, where each ``val``
            is a list [var1, var2, ...] of instance attribute names and/or method call
            arguments (of the decorated method). Values of the instance
            attributes/method call arguments will be saved to the location specified by
            ``key``.
    """
    property_names = property_names or {}

    def decorator(func):
        @wraps(func)
        def new_function(inst, *args, **kwargs):

            return_value = func(inst, *args, **kwargs)

            # Get a dict of the call arguments/keyword arguments by name:
            call_values = getcallargs(func, inst, *args, **kwargs)

            all_property_names = set()
            for names in property_names.values():
                all_property_names.update(names)

            property_values = {}
            for name in all_property_names:
                # If there is an instance attribute with that name, use that, otherwise
                # use the call value:
                if hasattr(inst, name):
                    property_values[name] = getattr(inst, name)
                else:
                    property_values[name] = call_values[name]

            # Save them:
            inst.set_properties(property_values, property_names)

            return return_value

        return new_function
    
    return decorator


def fastflatten(inarray, dtype):
    """A faster way of flattening our arrays than pylab.flatten.

    pylab.flatten returns a generator which takes a lot of time and memory
    to convert into a numpy array via array(list(generator)).  The problem
    is that generators don't know how many values they'll return until
    they're done. This algorithm produces a numpy array directly by
    first calculating what the length will be. It is several orders of
    magnitude faster. Note that we can't use numpy.ndarray.flatten here
    since our inarray is really a list of 1D arrays of varying length
    and/or single values, not a N-dimenional block of homogeneous data
    like a numpy array.

    Args:
        inarray (list): List of 1-D arrays to flatten.
        dtype (data-type): Type of the data in the arrays.

    Returns:
        :obj:`numpy:numpy.ndarray`: Flattened array.
    """
    total_points = np.sum([len(element) if np.iterable(element) else 1 for element in inarray])
    flat = np.empty(total_points,dtype=dtype)
    i = 0
    for val in inarray:
        if np.iterable(val):
            flat[i:i+len(val)] = val[:]
            i += len(val)
        else:
            flat[i] = val
            i += 1
    return flat


def max_or_zero(*args, **kwargs):
    """Returns max of the arguments or zero if sequence is empty.
    
    This protects the call to `max()` which would normally throw an error on an empty
    sequence.

    Args:
        *args: Items to compare.
        **kwargs: Passed to `max()`.

    Returns:
        : Max of *args.
    """
    if not args:
        return 0
    if not args[0]:
        return 0
    else:
        return max(*args, **kwargs)


def bitfield(arrays,dtype):
    """Converts a list of arrays of ones and zeros into a single
    array of unsigned ints of the given datatype.

    Args:
        arrays (list): List of numpy arrays consisting of ones and zeros.
        dtype (data-type): Type to convert to.

    Returns:
        :obj:`numpy:numpy.ndarray`: Numpy array with data type `dtype`.
    """
    n = {np.uint8: 8, np.uint16: 16, np.uint32: 32}
    if np.array_equal(arrays[0], 0):
        y = np.zeros(
            max([len(arr) if np.iterable(arr) else 1 for arr in arrays]), dtype=np.dtype
        )
    else:
        y = np.array(arrays[0], dtype=dtype)
    for i in range(1, n[dtype]):
        if np.iterable(arrays[i]):
            y |= arrays[i] << i
    return y


@contextlib.contextmanager
def suppress_mild_warnings(state=True):
    """A context manager which modifies compiler.suppress_mild_warnings

    Allows the user to suppress (or show) mild warnings for specific lines. Useful when
    you want to hide/show all warnings from specific lines.

    Arguments:
        state (bool): The new state for ``compiler.suppress_mild_warnings``. Defaults to
        ``True`` if not explicitly provided.
    """
    previous_warning_setting = compiler.suppress_mild_warnings
    compiler.suppress_mild_warnings = state
    yield
    compiler.suppress_mild_warnings = previous_warning_setting


@contextlib.contextmanager
def suppress_all_warnings(state=True):
    """A context manager which modifies compiler.suppress_all_warnings

    Allows the user to suppress (or show) all warnings for specific lines. Useful when
    you want to hide/show all warnings from specific lines.
    
    Arguments:
        state (bool): The new state for ``compiler.suppress_all_warnings``. Defaults to
        ``True`` if not explicitly provided.
    """
    previous_warning_setting = compiler.suppress_all_warnings
    compiler.suppress_all_warnings = state
    yield
    compiler.suppress_all_warnings = previous_warning_setting


class LabscriptError(Exception):
    """A *labscript* error.

    This is used to denote an error within the labscript suite itself.
    Is a thin wrapper of :obj:`Exception`.
    """
