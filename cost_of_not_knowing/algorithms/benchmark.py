# Author: Martin Schmidt
#         martin.schmidt@uni-trier.de
#         Mathias Sirvent
#         sirvent@posteo.de
import json
import logging
import time
from enum import Enum

import pyomo.environ as pyo
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition

from coding.cost_of_not_knowing.algorithms.algorithm import Algorithm
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.connections import CompressorStation
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.nodes import Entry, Exit
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.unit_helper import unit
from coding.cost_of_not_knowing.network.network_container import NetworkContainer
from coding.cost_of_not_knowing.network.pipe_container import PipeContainer
from coding.cost_of_not_knowing.utils.utils import (
    flow_unit,
    pressure_unit,
    norm_vol_per_hour_to_mass_flow,
    mass_flow_to_norm_vol_per_hour,
)

logger = logging.getLogger(__name__.split(".")[-1])


class Solver(Enum):
    ipopt = 1
    scip = 2


class Benchmark(Algorithm):
    """
    The benchmark algorithm on an MINLP.

    After initializing the class, the user can start the algorithm with the execute method.
    """

    def __init__(self, index: int, net_file: str, scn_file: str, solver: Solver):
        """
        Constructor of the algorithms model class.
        """
        start = time.time()
        super().__init__(index, net_file)
        self._nw_container = NetworkContainer(net_file, scn_file)
        self._model = None
        self._solver = solver
        end = time.time()
        self._read_network_runtime = end - start

    def execute(self):
        """
        Entry point to execute the benchmark algorithm.
        """
        start = time.time()
        self._build_model()
        obj, solver_runtime = self._solve_model()
        end = time.time()
        algorithm_runtime = end - start
        self._extract_solution_data_and_write(obj, algorithm_runtime, solver_runtime)

        # lp = self._result_directory.joinpath("lp.lp")
        # self._model.write(str(lp), io_options={"symbolic_solver_labels": True})

    def _build_model(self):
        """
        Setting up the network.
        """
        self._model = pyo.ConcreteModel()

        self._build_sets()
        self._build_variables()
        self._build_constraints()
        self._build_objective()

        logger.info("Building the pyomo network ... done")

    def _build_sets(self):
        """
        Building all sets required for the pyomo network.
        """
        self._model.arcs = [(u, v) for u, v in self._nw_container.network.edges()]
        self._model.pipes = [
            (u, v)
            for u, v in self._nw_container.network.edges()
            if isinstance(self._nw_container.network[u][v][0]["data"], PipeContainer)
        ]
        self._model.compressor_stations = [
            (u, v)
            for u, v in self._nw_container.network.edges()
            if isinstance(self._nw_container.network[u][v][0]["data"], CompressorStation)
        ]
        self._model.nodes = list(self._nw_container.network.nodes())
        self._model.bnd_nodes = [
            nodeId
            for nodeId in self._nw_container.network.nodes()
            if (
                isinstance(self._nw_container.network.nodes[nodeId]["data"], Entry)
                or isinstance(self._nw_container.network.nodes[nodeId]["data"], Exit)
            )
        ]

    def _build_variables(self):
        """
        Building the variables of the pyomo network. This method also takes care of unit conversions,
        so it's a dangerous place.
        """

        # Mass flow variables for arcs
        def arc_flow_bounds(m, u, v):
            """Rule to set the bounds for the flow variables at arcs."""
            lb_in_norm_vol_per_hour = self._nw_container.network[u][v][0]["data"].flow_min.to(flow_unit).magnitude
            ub_in_norm_vol_per_hour = self._nw_container.network[u][v][0]["data"].flow_max.to(flow_unit).magnitude
            lower_bound = norm_vol_per_hour_to_mass_flow(lb_in_norm_vol_per_hour)
            upper_bound = norm_vol_per_hour_to_mass_flow(ub_in_norm_vol_per_hour)
            return lower_bound, upper_bound

        self._model.arc_flow_vars = pyo.Var(self._model.arcs, bounds=arc_flow_bounds)

        # Pressure variables for nodes
        def pressure_bounds(m, node):
            """Rule to set the bounds for the pressure variables at nodes."""
            lower_bound = self._nw_container.network.nodes[node]["data"].pressure_min.to(pressure_unit).magnitude
            upper_bound = self._nw_container.network.nodes[node]["data"].pressure_max.to(pressure_unit).magnitude
            return lower_bound, upper_bound

        self._model.pressure_variables = pyo.Var(self._model.nodes, bounds=pressure_bounds)

        # Mass flow variables for boundary nodes
        def bnd_node_flow_bounds(m, bnd_node):
            """Rule to set the bounds for the flow variables at boundary nodes."""
            lb_in_norm_vol_per_hour = (
                self._nw_container.network.nodes[bnd_node]["data"].nomination_min.to(flow_unit).magnitude
            )
            ub_in_norm_vol_per_hour = (
                self._nw_container.network.nodes[bnd_node]["data"].nomination_max.to(flow_unit).magnitude
            )
            lower_bound = norm_vol_per_hour_to_mass_flow(lb_in_norm_vol_per_hour)
            upper_bound = norm_vol_per_hour_to_mass_flow(ub_in_norm_vol_per_hour)
            return lower_bound, upper_bound

        self._model.bnd_node_flow_variables = pyo.Var(self._model.bnd_nodes, bounds=bnd_node_flow_bounds)

        # Pressure increase variables for compressor stations
        def press_inc_bounds(m, u, v):
            """Rule to set the bounds for the flow variables at boundary nodes."""
            # bounds are interpreted in bar
            lower_bound = 5
            upper_bound = 30
            return lower_bound, upper_bound

        self._model.press_inc_variables = pyo.Var(self._model.compressor_stations, bounds=press_inc_bounds)

        logger.info("Building the pyomo variables ... done")

    def _build_constraints(self):

        # Flow balance
        def flow_balance_rule(m, node):
            """Rule to add flow balance constraints."""
            # node_data = self._nw_container.network.nodes[node]["data"]

            # out flow
            out_edges = self._nw_container.network.out_edges(node)
            out_flow = 0
            for (u, v) in out_edges:
                assert u == node
                out_flow += m.arc_flow_vars[(u, v)]

            # in flow
            in_edges = self._nw_container.network.in_edges(node)
            in_flow = 0
            for (u, v) in in_edges:
                assert v == node
                in_flow += m.arc_flow_vars[(u, v)]

            # rhs
            rhs = 0
            if isinstance(self._nw_container.network.nodes[node]["data"], Entry):
                rhs += m.bnd_node_flow_variables[node]
            elif isinstance(self._nw_container.network.nodes[node]["data"], Exit):
                rhs -= m.bnd_node_flow_variables[node]
            else:
                rhs = 0

            return out_flow - in_flow == rhs

        self._model.flow_balance_cons = pyo.Constraint(self._model.nodes, rule=flow_balance_rule)

        # Pressure loss
        def pressure_loss_rule(m, u, v):
            """
            Rule to add the pressure loss constraint for pipes.
            This method assumes pressure variables in bar and
            flow variables in mass flow (kg/s).
            """
            # data
            rhs = (
                self._nw_container.network[u][v][0]["data"].coefficient
                * m.arc_flow_vars[u, v]
                * abs(m.arc_flow_vars[u, v])
            )
            return m.pressure_variables[u] ** 2 - m.pressure_variables[v] ** 2 == rhs

        self._model.pressure_loss_cons = pyo.Constraint(self._model.pipes, rule=pressure_loss_rule)

        # Pressure increase at compressor stations
        def pressure_increase_rule(m, u, v):
            """
            Rule to add the pressure increase constraint for compressor stations.
            """
            return m.pressure_variables[v] == m.pressure_variables[u] + m.press_inc_variables[u, v]

        self._model.pressure_inc_cons = pyo.Constraint(self._model.compressor_stations, rule=pressure_increase_rule)

        logger.info("Building the pyomo constraints ... done")

    def _build_objective(self):
        """
        Building the objective function.
        """

        def objective_rule(m):
            """Objective rule"""
            return pyo.summation(self._model.press_inc_variables)

        self._model.objective = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

        logger.info("Building the pyomo objective ... done")

    def _solve_model(self) -> tuple:
        """Solving the network and return objective value and runtime."""

        fac = "scipampl" if self._solver == Solver.scip else "ipopt"
        solver = SolverFactory(fac)
        results = solver.solve(self._model, tee=True)  # noqa: F841
        if (results.solver.status == SolverStatus.ok) and (
            results.solver.termination_condition == TerminationCondition.optimal
        ):
            pass
        elif results.solver.termination_condition == TerminationCondition.infeasible:
            raise ValueError("SCIP returned status infeasible")
        else:
            raise ValueError(f"An unexpected error occured. SCIP returned status {results.solver.status}")
        return self._model.objective(), results.solver.Time

    def _extract_solution_data_and_write(self, obj: float, algorithm_runtime: float, solver_runtime: float):
        """Extracting data from the solution and write."""
        node_to_pressure = {}
        for node in self._model.nodes:
            node_to_pressure[node] = unit.Quantity(pyo.value(self._model.pressure_variables[node]), "bar")

        arcs_to_flows = {}
        for arc in self._model.arcs:
            norm_vol = mass_flow_to_norm_vol_per_hour(pyo.value(self._model.arc_flow_vars[arc]))
            arcs_to_flows[arc] = unit.Quantity(norm_vol, flow_unit)

        bnd_nodes_to_flows = {}
        for bnd_node in self._model.bnd_nodes:
            norm_vol = mass_flow_to_norm_vol_per_hour(pyo.value(self._model.bnd_node_flow_variables[bnd_node]))
            bnd_nodes_to_flows[bnd_node] = unit.Quantity(norm_vol, flow_unit)

        compressors_to_press_inc = {}
        for compressor_station in self._model.compressor_stations:
            compressors_to_press_inc[compressor_station] = unit.Quantity(
                pyo.value(self._model.press_inc_variables[compressor_station]),
                "bar",
            )

        self._write_solution(node_to_pressure, arcs_to_flows, bnd_nodes_to_flows, compressors_to_press_inc)
        dic = {
            "read_network_runtime": self._read_network_runtime,
            "algorithm_runtime": algorithm_runtime,
            "solver_runtime_inside_of_algorithm_runtime": solver_runtime,
            "total_runtime": self._read_network_runtime + algorithm_runtime,
            "objective": obj,
        }
        with open(str(self._result_directory.joinpath("summary.json")), "w") as f:
            f.write(json.dumps(dic, indent=4))
