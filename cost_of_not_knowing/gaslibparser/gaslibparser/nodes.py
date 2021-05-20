# Global imports
from collections import namedtuple

import numpy as np

# local imports
from .unit_helper import unit

_node_defaults = (  # node_id is required
    None,  # pos
    0,  # height
    unit.Quantity(0, "bar"),  # pressure_min
    unit.Quantity(np.Infinity, "bar"),  # pressure_max
)


class Node(
    namedtuple("NodeNamedTuple", ["node_id", "pos", "height", "pressure_min", "pressure_max"], defaults=_node_defaults),
):
    """
    Node in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


_entry_defaults = (
    _node_defaults
    + (
        unit.Quantity(0, "m_cube_per_hour"),  # flow_min
        unit.Quantity(np.Infinity, "m_cube_per_hour"),  # flow_max
        unit.Quantity(0, "m_cube_per_hour"),  # nomination_min
        unit.Quantity(0, "m_cube_per_hour"),  # nomination_max
    )
    + (None,) * 9
)


class Entry(
    namedtuple(
        "EntryNamedTuple",
        Node._fields
        + (
            "flow_min",
            "flow_max",
            "nomination_min",
            "nomination_max",
            "gas_temp",
            "calorific_value",
            "norm_density",
            "heat_coeff_A",
            "heat_coeff_B",
            "heat_coeff_C",
            "molar_mass",
            "pseudocritical_pressure",
            "pseudocritical_temperature",
        ),
        defaults=_entry_defaults,
    ),
):
    """
    Entry node in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """


_exit_defaults = _node_defaults + (
    unit.Quantity(0, "m_cube_per_hour"),  # flow_min
    unit.Quantity(np.Infinity, "m_cube_per_hour"),  # flow_max
    unit.Quantity(0, "m_cube_per_hour"),  # nomination_min
    unit.Quantity(0, "m_cube_per_hour"),  # nomination_max
)


class Exit(
    namedtuple(
        "ExitNamedTuple",
        Node._fields + ("flow_min", "flow_max", "nomination_min", "nomination_max"),
        defaults=_exit_defaults,
    ),
):
    """
    Exit node in a gas transport network.
    Realized by inheritance from a suitable namedtuple.
    """
