from __future__ import annotations

import math

import numpy as np


class SamplingPoint:
    def __init__(
        self,
        iteration: int,
        pressure_in: float,  # bar
        pressure_out: float,  # bar
        mass_flow: float,  # kg_per_second,
        pressure_in_fixed: bool,
        pressure_out_fixed: bool,
        flow_fixed: bool,
        coefficient: float,
        lipschitz: float,
    ):
        self.iteration = iteration

        self.pressure_in = pressure_in
        self.pressure_out = pressure_out
        self.mass_flow = mass_flow

        self.pressure_in_fixed = pressure_in_fixed
        self.pressure_out_fixed = pressure_out_fixed
        self.flow_fixed = flow_fixed

        self.F = None  # dimensionless
        self.pressure_in_relative_error = None  # dimensionless
        self.pressure_out_relative_error = None  # dimensionless
        self.mass_flow_relative_error = None  # dimensionless
        self.pressure_in_lb = None  # float in bar
        self.pressure_in_ub = None  # float in bar
        self.pressure_out_lb = None  # float in bar
        self.pressure_out_ub = None  # float in bar
        self.mass_flow_lb = None  # float in kg_per_second
        self.mass_flow_ub = None  # float in kg_per_second

        self._set_f_value(coefficient)
        self._set_relative_errors(coefficient)
        self._set_bounds(lipschitz)

    def is_excluded(self, other: SamplingPoint) -> bool:
        """
        Checks if the object is excluded by other SamplingPoint.
        """
        if (
            self.pressure_in_fixed != other.pressure_in_fixed
            or self.pressure_out_fixed != other.pressure_out_fixed
            or self.flow_fixed != other.flow_fixed
        ):
            assert "Sampling points cannot be compared"
        if not self.pressure_in_fixed and (
            self.pressure_in < other.pressure_in_lb or self.pressure_in > other.pressure_in_ub
        ):
            return False
        if not self.pressure_out_fixed and (
            self.pressure_out < other.pressure_out_lb or self.pressure_out > other.pressure_out_ub
        ):
            return False
        if not self.flow_fixed and (self.mass_flow < other.mass_flow_lb or self.mass_flow > other.mass_flow_ub):
            return False
        return True

    def _set_f_value(self, coefficient: float):
        """
        Calculates and sets the F value for a sampling point. F is dimensionless.
        """
        lhs = self.pressure_in ** 2 - self.pressure_out ** 2
        rhs = coefficient * self.mass_flow * abs(self.mass_flow)
        self.F = lhs - rhs

    def _set_relative_errors(self, coefficient: float):
        """
        Calculates and sets the relative errors for a sampling point.
        """
        # F = self.pressure_in ** 2 - self.pressure_out ** 2 - coefficient * self.mass_flow * abs(self.mass_flow)
        pressure_in_squared = self.pressure_out ** 2 + coefficient * self.mass_flow * abs(self.mass_flow)
        if pressure_in_squared >= 0:
            pressure_in = math.sqrt(pressure_in_squared)
            self.pressure_in_relative_error = abs(self.pressure_in - pressure_in) / self.pressure_in
        else:
            self.pressure_in_relative_error = np.inf

        pressure_out_squared = self.pressure_in ** 2 - coefficient * self.mass_flow * abs(self.mass_flow)
        if pressure_out_squared >= 0:
            pressure_out = math.sqrt(pressure_out_squared)
            self.pressure_out_relative_error = abs(self.pressure_out - pressure_out) / self.pressure_out
        else:
            self.pressure_out_relative_error = np.inf

        mass_flow_squared = abs((self.pressure_in ** 2 - self.pressure_out ** 2) / coefficient)
        if mass_flow_squared >= 0 and self.mass_flow != 0.0:
            mass_flow = math.sqrt(mass_flow_squared)
            self.mass_flow_relative_error = abs(abs(self.mass_flow) - mass_flow) / abs(self.mass_flow)
        else:
            self.mass_flow_relative_error = np.inf

    def _set_bounds(self, lipschitz: float):
        """
        Set the lower bounds and upper bounds for the boxes to exclude.
        """
        box_size = abs(self.F) / lipschitz
        self.pressure_in_lb = self.pressure_in - box_size if not self.pressure_in_fixed else self.pressure_in_lb
        self.pressure_in_ub = self.pressure_in + box_size if not self.pressure_in_fixed else self.pressure_in_ub
        self.pressure_out_lb = self.pressure_out - box_size if not self.pressure_out_fixed else self.pressure_out_lb
        self.pressure_out_ub = self.pressure_out + box_size if not self.pressure_out_fixed else self.pressure_out_ub
        self.mass_flow_lb = self.mass_flow - box_size if not self.flow_fixed else self.mass_flow_lb
        self.mass_flow_ub = self.mass_flow + box_size if not self.flow_fixed else self.mass_flow_ub
