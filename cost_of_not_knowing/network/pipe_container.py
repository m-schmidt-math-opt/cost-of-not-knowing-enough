# Author: Martin Schmidt
#         martin.schmidt@uni-trier.de
#         Mathias Sirvent
#         sirvent@posteo.de
from __future__ import annotations

import itertools
import math
import time
from typing import Union, Optional

import numpy as np

from coding.cost_of_not_knowing.gaslibparser.gaslibparser.connections import Pipe
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.nodes import Entry, Node, Exit
from coding.cost_of_not_knowing.network.sampling_point import SamplingPoint
from coding.cost_of_not_knowing.utils.utils import (
    length_unit,
    nikuradse,
    speed_of_sound,
    norm_vol_per_hour_to_mass_flow,
    pressure_unit,
    flow_unit,
)


class PipeContainer:
    def __init__(
        self,
        pipe: Pipe,
        from_node_data: Union[Entry, Node, Exit],
        to_node_data: Union[Entry, Node, Exit],
        sampling_steps: Optional[tuple],
        is_box_algorithm: bool,
    ):
        self.arc_id = pipe.arc_id
        self.from_node = pipe.from_node
        self.to_node = pipe.to_node
        self.flow_min = pipe.flow_min  # Quantity in m_cube_per_hour
        self.flow_max = pipe.flow_max  # Quantity in m_cube_per_hour
        self.length = pipe.length
        self.diameter = pipe.diameter
        self.roughness = pipe.roughness
        self.pressure_max = pipe.pressure_max
        self.heat_transfer_coefficient = pipe.heat_transfer_coefficient

        self.coefficient = None
        self._set_coefficient()

        self.lipschitz = None
        self.pressure_in_fixed = None  # boolean if pressure_in is already fixed because of tight bounds
        self.pressure_out_fixed = None  # boolean if pressure_out is already fixed because of tight bounds
        self.flow_fixed = None  # boolean if flow is already fixed because of tight bounds
        self.successful_points = []  # feasible points
        self.sampling_points_pending = []  # sampling points that have NOT been added to the optimization problem
        self.sampling_points_build_in = []  # sampling points that have been added to the optimization problem

        is_box_algorithm and self._set_variable_fixes(from_node_data, to_node_data)
        is_box_algorithm and self._set_lipschitz(from_node_data, to_node_data)
        start = time.time()
        is_box_algorithm and self._add_initial_sampling_points(from_node_data, to_node_data, sampling_steps)
        end = time.time()
        self.sampling_runtime = end - start

    def get_f(
        self,
        pressure_in: float,  # bar,
        pressure_out: float,  # bar,
        mass_flow: float,  # kg_per_second,
    ) -> float:
        """
        Returns the F value for a solution triple.
        """
        sp_helper = SamplingPoint(
            0,
            pressure_in,
            pressure_out,
            mass_flow,
            self.pressure_in_fixed,
            self.pressure_out_fixed,
            self.flow_fixed,
            self.coefficient,
            self.lipschitz,
        )
        return sp_helper.F

    def get_relative_error(
        self,
        pressure_in: float,  # bar,
        pressure_out: float,  # bar,
        mass_flow: float,  # kg_per_second,
    ) -> float:
        """
        Returns the relative error for a solution triple.
        """
        sp_helper = SamplingPoint(
            0,
            pressure_in,
            pressure_out,
            mass_flow,
            self.pressure_in_fixed,
            self.pressure_out_fixed,
            self.flow_fixed,
            self.coefficient,
            self.lipschitz,
        )
        return max(
            sp_helper.pressure_in_relative_error,
            sp_helper.pressure_out_relative_error,
            sp_helper.mass_flow_relative_error,
        )

    def add_successful_point(
        self,
        iteration: int,
        pressure_in: float,  # bar,
        pressure_out: float,  # bar,
        mass_flow: float,  # kg_per_second,
    ):
        """
        Add a successful point.
        """
        sp = SamplingPoint(
            iteration,
            pressure_in,
            pressure_out,
            mass_flow,
            self.pressure_in_fixed,
            self.pressure_out_fixed,
            self.flow_fixed,
            self.coefficient,
            self.lipschitz,
        )
        self.successful_points.append(sp)

    def add_sampling_point(
        self,
        iteration: int,
        pressure_in: float,  # bar,
        pressure_out: float,  # bar,
        mass_flow: float,  # kg_per_second,
        add_without_check: bool = False,
    ):
        """
        Add a sampling point.
        """
        sp = SamplingPoint(
            iteration,
            pressure_in,
            pressure_out,
            mass_flow,
            self.pressure_in_fixed,
            self.pressure_out_fixed,
            self.flow_fixed,
            self.coefficient,
            self.lipschitz,
        )
        if add_without_check:
            self.sampling_points_pending.append(sp)
        else:
            sps = self.sampling_points_pending + [sp]

            # Sort sampling points by F
            sps.sort(key=lambda x: abs(x.F), reverse=True)

            # Check if the sampling point has already been removed by a sampling point with a bigger box.
            # In this case we do not need the sampling point.
            index = sps.index(sp)
            for sp_index in range(index):
                if sp.is_excluded(sps[sp_index]):
                    return

            # Check if the new sampling point removes sampling points with smaller boxes. In this case, remove
            # these sampling points and update the sampling_points_pending.
            sps_copy = sps[:]
            indexes_to_exclude = []
            for sp_index in range(index + 1, len(sps)):
                if sps_copy[sp_index].is_excluded(sp):
                    indexes_to_exclude.append(sp_index)
            for ele in sorted(indexes_to_exclude, reverse=True):
                del sps[ele]
            self.sampling_points_pending = sps

    def _set_coefficient(self):
        """
        Set the pipe coefficient. Holds for pressure in bar and mass flow in kg/s.
        """
        diameter = self.diameter.to(length_unit).magnitude
        length = self.length.to(length_unit).magnitude
        roughness = self.roughness.to(length_unit).magnitude
        friction_coeff = nikuradse(diameter, roughness)
        self.coefficient = 1e-10 * (16 * friction_coeff * speed_of_sound ** 2 * length) / (diameter ** 5 * math.pi ** 2)

    def _set_lipschitz(
        self,
        from_node_data: Union[Entry, Node, Exit],
        to_node_data: Union[Entry, Node, Exit],
    ):
        """
        Set global lipschitz constant. Note that lipschitz constant is dimensionless and depends on the selected units.
        """
        p_in_max = from_node_data.pressure_max.to(pressure_unit).magnitude
        p_out_max = to_node_data.pressure_max.to(pressure_unit).magnitude
        mass_flow_min = norm_vol_per_hour_to_mass_flow(self.flow_min.to(flow_unit).magnitude)
        mass_flow_max = norm_vol_per_hour_to_mass_flow(self.flow_max.to(flow_unit).magnitude)
        q_max = max(abs(mass_flow_min), abs(mass_flow_max))
        if self.flow_fixed:
            self.lipschitz = 2.0 * p_in_max + 2.0 * p_out_max
        else:
            self.lipschitz = 2.0 * p_in_max + 2.0 * p_out_max + 2.0 * self.coefficient * q_max

    def _set_variable_fixes(
        self,
        from_node_data: Union[Entry, Node, Exit],
        to_node_data: Union[Entry, Node, Exit],
    ):
        """
        Sets the variables to check if flow or pressure is already fixed due to tight bounds.
        """
        self.flow_fixed = math.isclose(
            norm_vol_per_hour_to_mass_flow(self.flow_min.to(flow_unit).magnitude),
            norm_vol_per_hour_to_mass_flow(self.flow_max.to(flow_unit).magnitude),
        )
        self.pressure_in_fixed = math.isclose(
            from_node_data.pressure_min.to(pressure_unit).magnitude,
            from_node_data.pressure_max.to(pressure_unit).magnitude,
        )
        self.pressure_out_fixed = math.isclose(
            to_node_data.pressure_min.to(pressure_unit).magnitude,
            to_node_data.pressure_max.to(pressure_unit).magnitude,
        )

    def _add_initial_sampling_points(
        self,
        from_node_data: Union[Entry, Node, Exit],
        to_node_data: Union[Entry, Node, Exit],
        sampling_steps: Optional[tuple],
    ):
        """
        Adding initial sampling points.
        """
        if not sampling_steps:
            return
        pressure_steps = sampling_steps[0]  # in bar
        mass_flow_steps = sampling_steps[1]  # in kg_per_second

        pressure_in_min = from_node_data.pressure_min.to(pressure_unit).magnitude
        pressure_in_max = from_node_data.pressure_max.to(pressure_unit).magnitude
        pressure_out_min = to_node_data.pressure_min.to(pressure_unit).magnitude
        pressure_out_max = to_node_data.pressure_max.to(pressure_unit).magnitude
        mass_flow_min = norm_vol_per_hour_to_mass_flow(self.flow_min.to(flow_unit).magnitude)
        mass_flow_max = norm_vol_per_hour_to_mass_flow(self.flow_max.to(flow_unit).magnitude)

        pressure_ins = np.arange(pressure_in_min, pressure_in_max, pressure_steps).tolist()
        pressure_ins = pressure_ins if pressure_ins else [(pressure_in_min + pressure_in_max) * 0.5]
        pressure_outs = np.arange(pressure_out_min, pressure_out_max, pressure_steps).tolist()
        pressure_outs = pressure_outs if pressure_outs else [(pressure_out_min + pressure_out_max) * 0.5]
        mass_flows = np.arange(mass_flow_min, mass_flow_max, mass_flow_steps).tolist()
        mass_flows = mass_flows if mass_flows else [(mass_flow_min + mass_flow_max) * 0.5]

        [self.add_sampling_point(-1, *sp) for sp in itertools.product(pressure_ins, pressure_outs, mass_flows)]
