Connection Table
================
The connection table maps out the way input/output devices are connected to each other in your lab, and the channels (individual inputs/outputs) they have. The devices in your lab should be connected in a similar way to that shown in the figure below.

TODO: insert figure!

Here we see two :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>` instances in the top tier of the diagram. They do not have a parent device that tells them when to update their output (this is true for all :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>` instances). However, all but one (the master pseudoclock device) must be triggered by an output clocked by the master pseudoclock device. 

Each :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>` instance should have one or more :py:class:`Pseudoclock <labscript.Pseudoclock>` children. Some :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>` instances may automatically create these children for you (check the device specific documentation). In turn, each :py:class:`Pseudoclock <labscript.Pseudoclock>` will have one of more :py:class:`ClockLine <labscript.ClockLine>` instances connected to it. These :py:class:`ClockLine <labscript.ClockLine>` instances generally refer to physical outputs of a device which will be used to clock another device. However, in some cases, one or more :py:class:`ClockLine <labscript.ClockLine>` instances may be internally created for you (check the device specific documentation).

If a device is not a :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>`, it must be connected to one via a clockline. such devices inherit from :py:class:`IntermediateDevice <labscript.IntermediateDevice>`. Inputs and outputs are then connected to these devices. If a :py:class:`PseudoclockDevice <labscript.PseudoclockDevice>` also has outputs that are not used for a :py:class:`ClockLine <labscript.ClockLine>`, then an :py:class:`IntermediateDevice <labscript.IntermediateDevice>` is internally instantiated, and should be made available through the ``PseudoclockDevice.direct_outputs`` attribute (for example see PulseBlaster implementation TODO: link!).

