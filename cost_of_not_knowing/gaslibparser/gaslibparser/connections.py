# Global imports
from collections import namedtuple


class Arc(namedtuple("ArcNamedTuple", ["arc_id", "from_node", "to_node", "flow_min", "flow_max"])):
    """
    Abstract arc in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class Pipe(
    namedtuple(
        "PipeNamedTuple",
        Arc._fields
        + (
            "length",
            "diameter",
            "roughness",
            "pressure_max",
            "heat_transfer_coefficient",
        ),
        defaults=(None,) * 2,
    ),
):
    """
    Pipe of a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class CompressorStation(
    namedtuple(
        "CompressorStationNamedTuple",
        Arc._fields
        + (
            "pressure_differential_min",
            "pressure_differential_max",
            "flow_threshold",
            "pressure_in_min",
            "pressure_out_max",
            "max_pressure_increase",
            "max_pressure_ratio",
            "min_pressure_increase",
            "min_pressure_ratio",
        ),
        defaults=(None,) * 9,
    ),
):
    """
    Compressor station of a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class ControlValve(
    namedtuple(
        "ControlValveNamedTuple",
        Arc._fields
        + (
            "pressure_differential_min",
            "pressure_differential_max",
            "flow_threshold",
            "pressure_in_min",
            "pressure_out_max",
        ),
        defaults=(None,) * 5,
    ),
):
    """
    Control valve of a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class Valve(namedtuple("ValveNamedTuple", Arc._fields + ("pressure_differential_max",), defaults=(None,) * 1)):
    """
    Valve of a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class ShortPipe(namedtuple("ShortPipeNamedTuple", Arc._fields)):
    """
    Short pipe in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


class Resistor(namedtuple("ResistorNamedTuple", Arc._fields + ("diameter", "drag_factor"), defaults=(None,) * 2)):
    """
    Resistor in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """
