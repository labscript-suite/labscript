"""Classes for device channels that are inputs"""

from .base import Device
from .utils import set_passed_properties


class AnalogIn(Device):
    """Analog Input for use with all devices that have an analog input."""
    description = "Analog Input"

    @set_passed_properties(property_names={})
    def __init__(
        self, name, parent_device, connection, scale_factor=1.0, units="Volts", **kwargs
    ):
        """Instantiates an Analog Input.

        Args:
            name (str): python variable to assign this input to.
            parent_device (:obj:`IntermediateDevice`): Device input is connected to.
            scale_factor (float, optional): Factor to scale the recorded values by.
            units (str, optional): Units of the input.
            **kwargs: Keyword arguments passed to :func:`Device.__init__`.
        """
        self.acquisitions = []
        self.scale_factor = scale_factor
        self.units=units
        Device.__init__(self, name, parent_device, connection, **kwargs)

    def acquire(
        self, label, start_time, end_time, wait_label="", scale_factor=None, units=None
    ):
        """Command an acquisition for this input.

        Args:
            label (str): Unique label for the acquisition. Used to identify the saved trace.
            start_time (float): Time, in seconds, when the acquisition should start.
            end_time (float): Time, in seconds, when the acquisition should end.
            wait_label (str, optional): 
            scale_factor (float): Factor to scale the saved values by.
            units: Units of the input, consistent with the unit conversion class.

        Returns:
            float: Duration of the acquistion, equivalent to `end_time - start_time`.
        """
        if scale_factor is None:
            scale_factor = self.scale_factor
        if units is None:
            units = self.units
        self.acquisitions.append(
            {
                "start_time": start_time,
                "end_time": end_time,
                "label": label,
                "wait_label": wait_label,
                "scale_factor": scale_factor,
                "units": units,
            }
        )
        return end_time - start_time

