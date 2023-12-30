from inspect import getcallargs
from functools import wraps

_RemoteConnection = None
ClockLine = None
PseudoClockDevice = None


def is_remote_connection(connection):
    """Returns whether the connection is an instance of ``_RemoteConnection``
    
    This function defers and caches the import of ``_RemoteConnection``. This both
    breaks the circular import between ``Device`` and ``_RemoteConnection``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    if _RemoteConnection is None:
        from .remote import _RemoteConnection
    return isinstance(connection, _RemoteConnection)


def is_clock_line(device):
    """Returns whether the connection is an instance of ``ClockLine``
    
    This function defers and caches the import of ``ClockLine``. This both
    breaks the circular import between ``Device`` and ``ClockLine``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    if ClockLine is None:
        from .labscript import ClockLine
    return isinstance(device, _RemoteConnection)


def is_pseudoclock_device(device):
    """Returns whether the connection is an instance of ``PseudoclockDevice``
    
    This function defers and caches the import of ``_RemoteConnection``. This both
    breaks the circular import between ``Device`` and ``_RemoteConnection``, while
    maintaining reasonable performance (this performs better than importing each time as
    the lookup in the modules hash table is slower).
    """
    if PseudoclockDevice is None:
        from .labscript import PseudoclockDevice
    return isinstance(device, PseudoclockDevice)


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


class LabscriptError(Exception):
    """A *labscript* error.

    This is used to denote an error within the labscript suite itself.
    Is a thin wrapper of :obj:`Exception`.
    """
