# Author: Martin Schmidt
#         martin.schmidt@uni-trier.de
#         Mathias Sirvent
#         sirvent@posteo.de

from __future__ import annotations

import json
import logging
import time
from typing import List, Optional

import gurobipy as gp
import matplotlib.pyplot as plt
import numpy as np
from gurobipy import LinExpr, GRB
from matplotlib.patches import Rectangle

from coding.cost_of_not_knowing.algorithms.algorithm import Algorithm
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.connections import CompressorStation
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.nodes import Exit, Entry
from coding.cost_of_not_knowing.gaslibparser.gaslibparser.unit_helper import unit
from coding.cost_of_not_knowing.network.network_container import NetworkContainer
from coding.cost_of_not_knowing.network.pipe_container import PipeContainer
from coding.cost_of_not_knowing.utils.utils import (
    norm_vol_per_hour_to_mass_flow,
    flow_unit,
    pressure_unit,
    mass_flow_to_norm_vol_per_hour,
)

logger = logging.getLogger(__name__.split(".")[-1])


class BoxAlgorithm(Algorithm):
    """
    The box algorithm.

    After initializing the class, the user can start the algorithm with the execute method.

    """

    def __init__(
        self,
        index: int,
        net_file: str,
        scn_file: str,
        sampling_steps: Optional[tuple],  # first: bar steps for pressure, second: mass flow steps in kg_per_second
        error_value: float,
        absolute_error: bool,
    ):
        """
        Constructor of the box algorithm class.
        """
        start = time.time()
        super().__init__(index, net_file)
        self._nw_container = NetworkContainer(net_file, scn_file, sampling_steps, True)
        self._sampling_steps = sampling_steps
        self._error_value = error_value
        self._absolute_error = absolute_error
        self._model = gp.Model("box_algorithm")
        self._model.setParam("LogToConsole", 0)
        self._model.setParam("LogFile", (str(self._result_directory.joinpath("gurobi.log"))))
        self._sol_dir = self._result_directory.joinpath("solutions")
        self._sol_dir.mkdir(parents=True, exist_ok=True)
        self._lp_dir = self._result_directory.joinpath("lps")
        self._lp_dir.mkdir(parents=True, exist_ok=True)
        self._pipe_dir = self._result_directory.joinpath("pipe_stats")
        self._pipe_dir.mkdir(parents=True, exist_ok=True)
        end = time.time()
        self._read_network_runtime = end - start

    def execute(self):
        """
        Entry point to execute the box algorithm.
        """
        logger.info("--------- Start Iteration 0 only with initial Sampling Points ---------")
        start = time.time()
        solver_runtime = 0.0
        iteration = 0
        self._build_init_model()
        feasible, runtime, obj = self._solve(iteration)
        solver_runtime += runtime

        optimal, max_error = self._is_optimal(iteration)
        time_limit = False
        while feasible and not optimal:
            iteration += 1
            logger.info(f"--------- Start Iteration {iteration} ---------")
            self._build_sampling_point_constraints()
            feasible, runtime, obj = self._solve(iteration)
            solver_runtime += runtime
            optimal, max_error = self._is_optimal(iteration)
            elapsed_time = time.time() - start
            logger.info(f"--------- elapsed time {elapsed_time} seconds ---------")
            if elapsed_time >= 900.0:
                logger.info("Time limit reached")
                time_limit = True
                break
            logger.info(f"--------- missing {900 - elapsed_time} seconds until time limit ---------")
            logger.info("")
        end = time.time()
        algorithm_runtime = end - start

        self._write_config(time_limit)
        if not time_limit:
            self._write_out_pipe_data()
            self._write_box_solution()
            self._write_summary(iteration, obj, max_error, algorithm_runtime, solver_runtime)

    def _solve(self, iteration: int) -> tuple:
        """
        Wrapper to write out the LP file, to solve it, and to write out the solution. Returns True, if feasible and
        False if infeasible. Raise an ValueError otherwise. Moreover, it returns the runtime and the objective value.
        """
        self._model.write(str(self._lp_dir.joinpath(f"lp_{iteration}.lp")))
        self._model.optimize()
        runtime = self._write_solution_stats(iteration)
        if self._model.Status == GRB.OPTIMAL:
            logger.info("Feasible")
            return True, runtime, self._model.objVal
        elif self._model.Status == GRB.INFEASIBLE:
            logger.info("Infeasible")
            return False, runtime, np.inf
        raise ValueError(f"An unexpected error occured. Gurobi returned status {self._model.Status}")

    def _is_optimal(self, iteration: int) -> tuple:
        """
        Checks the optimality of the current solution. Moreover, the method adds new sampling points if optimality is
        not reached yet. Returns True, if optimality is reached and False otherwise. Moreover, it returns the maximum
        relative or absolute error depending on the configuration.
        """
        optimality = True
        errors = []
        for (u, v) in self._nw_container.network.edges():
            pipe = self._nw_container.network[u][v][0]["data"]
            if isinstance(pipe, PipeContainer):
                mass_flow_sol = self._edge_flow_vars[pipe.arc_id].x
                pressure_in_sol = self._node_press_vars[u].x
                pressure_out_sol = self._node_press_vars[v].x
                optimal, error = self._nw_container.handle_new_point(
                    iteration,
                    pipe,
                    pressure_in_sol,
                    pressure_out_sol,
                    mass_flow_sol,
                    self._error_value,
                    self._absolute_error,
                )
                errors.append(error)
                optimality = False if not optimal else optimality
        return optimality, max(errors)  # noqa: R504

    def _build_init_model(self):
        """
        Build init model.
        """

        self._build_data()
        self._build_variables()
        self._build_constraints()
        self._build_sampling_point_constraints()

    def _build_data(self):
        """
        Build data.
        """

        self._edge_flow_data = {}
        self._node_pressure_data = {}
        self._node_flow_data = {}
        self._compr_incr_data = {}

        for (u, v) in self._nw_container.network.edges():
            edge_data = self._nw_container.network[u][v][0]["data"]
            edge_id = edge_data.arc_id
            lb = norm_vol_per_hour_to_mass_flow(edge_data.flow_min.to(flow_unit).magnitude)
            ub = norm_vol_per_hour_to_mass_flow(edge_data.flow_max.to(flow_unit).magnitude)
            self._edge_flow_data.update({edge_id: [lb, ub]})
            if isinstance(edge_data, CompressorStation):
                lb = 5
                ub = 30
                self._compr_incr_data.update({edge_id: [lb, ub]})

        for nd in self._nw_container.network.nodes():
            node_data = self._nw_container.network.nodes[nd]["data"]
            node_id = node_data.node_id
            lb = node_data.pressure_min.to(pressure_unit).magnitude
            ub = node_data.pressure_max.to(pressure_unit).magnitude
            self._node_pressure_data.update({node_id: [lb, ub]})
            if isinstance(node_data, Entry) or isinstance(node_data, Exit):
                lb_in_norm_vol_per_hour = node_data.nomination_min.to(flow_unit).magnitude
                ub_in_norm_vol_per_hour = node_data.nomination_max.to(flow_unit).magnitude
                lb = norm_vol_per_hour_to_mass_flow(lb_in_norm_vol_per_hour)
                ub = norm_vol_per_hour_to_mass_flow(ub_in_norm_vol_per_hour)
                self._node_flow_data.update({node_id: [lb, ub]})

    def _build_variables(self):
        """
        Build variables.
        """

        edge_flow, lbs, ubs = gp.multidict(self._edge_flow_data)
        self._edge_flow_vars = self._model.addVars(edge_flow, name="flow", lb=lbs, ub=ubs)
        node_pressure, lbs, ubs = gp.multidict(self._node_pressure_data)
        self._node_press_vars = self._model.addVars(node_pressure, name="pressure", lb=lbs, ub=ubs)
        node_flow, lbs, ubs = gp.multidict(self._node_flow_data)
        self._node_flow_vars = self._model.addVars(node_flow, name="flow", lb=lbs, ub=ubs)
        if self._compr_incr_data:
            compressor_increase, lbs, ubs = gp.multidict(self._compr_incr_data)
            self._compr_vars = self._model.addVars(compressor_increase, name="compressor_incr", lb=lbs, ub=ubs, obj=1.0)
        else:
            self._compr_vars = {}

    def _build_constraints(self):
        """
        Build constraints.
        """

        # Flow conservation
        for node in self._nw_container.network.nodes():

            # out flow
            out_flow = LinExpr()
            for (u, v) in self._nw_container.network.out_edges(node):
                out_flow.add(self._edge_flow_vars[self._nw_container.network[u][v][0]["data"].arc_id])

            # in flow
            in_flow = LinExpr()
            for (u, v) in self._nw_container.network.in_edges(node):
                in_flow.add(self._edge_flow_vars[self._nw_container.network[u][v][0]["data"].arc_id])

            # rhs
            rhs = LinExpr()
            node_id = self._nw_container.network.nodes[node]["data"].node_id
            if isinstance(self._nw_container.network.nodes[node]["data"], Entry):
                rhs.add(self._node_flow_vars[node_id], 1.0)
            elif isinstance(self._nw_container.network.nodes[node]["data"], Exit):
                rhs.add(self._node_flow_vars[node_id], -1.0)
            else:
                rhs = 0

            self._model.addConstr(out_flow - in_flow == rhs, f"flow_conservation_{node_id}")

        # Pressure increase constraint for compressor stations
        for (u, v) in self._nw_container.network.edges():
            edge_data = self._nw_container.network[u][v][0]["data"]
            if isinstance(edge_data, CompressorStation):
                edge_id = edge_data.arc_id
                self._model.addConstr(
                    self._node_press_vars[v] == self._node_press_vars[u] + self._compr_vars[edge_id],
                    f"compressor_increase_{edge_id}",
                )

    def _build_sampling_point_constraints(self):
        """
        Build pending sampling point constraints.
        """
        for (u, v) in self._nw_container.network.edges():
            edge_data = self._nw_container.network[u][v][0]["data"]
            if isinstance(edge_data, PipeContainer):
                arc_id = edge_data.arc_id
                for sp in edge_data.sampling_points_pending:
                    z = []
                    it = sp.iteration
                    if not sp.pressure_in_fixed:
                        tag = "pressure_in"
                        node_lb = self._nw_container.network.nodes[u]["data"].pressure_min.to(pressure_unit).magnitude
                        node_ub = self._nw_container.network.nodes[u]["data"].pressure_max.to(pressure_unit).magnitude
                        sp_lb = sp.pressure_in_lb
                        sp_ub = sp.pressure_in_ub
                        pressure_in_var = self._node_press_vars[u]
                        z.extend(
                            self._add_box(id(sp), it, arc_id, tag, pressure_in_var, node_lb, node_ub, sp_lb, sp_ub),
                        )
                    if not sp.pressure_out_fixed:
                        tag = "pressure_out"
                        node_lb = self._nw_container.network.nodes[v]["data"].pressure_min.to(pressure_unit).magnitude
                        node_ub = self._nw_container.network.nodes[v]["data"].pressure_max.to(pressure_unit).magnitude
                        sp_lb = sp.pressure_out_lb
                        sp_ub = sp.pressure_out_ub
                        pressure_out_var = self._node_press_vars[v]
                        z.extend(
                            self._add_box(id(sp), it, arc_id, tag, pressure_out_var, node_lb, node_ub, sp_lb, sp_ub),
                        )
                    if not sp.flow_fixed:
                        tag = "mass_flow"
                        node_lb = self._nw_container.network[u][v][0]["data"].flow_min.to(flow_unit).magnitude
                        node_ub = self._nw_container.network[u][v][0]["data"].flow_max.to(flow_unit).magnitude
                        sp_lb = sp.mass_flow_lb
                        sp_ub = sp.mass_flow_ub
                        mass_flow_var = self._edge_flow_vars[arc_id]
                        z.extend(self._add_box(id(sp), it, arc_id, tag, mass_flow_var, node_lb, node_ub, sp_lb, sp_ub))
                    self._model.addConstr(
                        LinExpr([1.0 for _ in z], z) >= 1,
                        name=f"cover_{arc_id}_iteration_{it}_sampling_point_id_{id(sp)}",
                    )
                edge_data.sampling_points_build_in.extend(edge_data.sampling_points_pending)
                edge_data.sampling_points_pending = []

    def _add_box(
        self,
        unique_id: int,
        it: int,
        pipe_id: str,
        tag: str,
        var: gp.Var,
        lb: float,
        ub: float,
        sp_lb: float,
        sp_ub: float,
    ) -> List:
        """
        Add the box constraints.
        """
        z_1 = self._model.addVar(
            vtype=gp.GRB.BINARY,
            name=f"z_1_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        z_2 = self._model.addVar(
            vtype=gp.GRB.BINARY,
            name=f"z_2_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        self._model.addConstr(
            var - (ub - sp_ub) * z_1 <= sp_ub,
            name=f"con_1_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        self._model.addConstr(
            var - (sp_ub - lb) * z_1 >= lb,
            name=f"con_2_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        self._model.addConstr(
            var - (ub - sp_lb) * (1 - z_2) <= sp_lb,
            name=f"con_3_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        self._model.addConstr(
            var - (sp_lb - lb) * (1 - z_2) >= lb,
            name=f"con_4_{pipe_id}_{tag}_iter_{it}_sampling_point_id_{unique_id}",
        )
        return [z_1, z_2]

    def _write_out_pipe_data(self):
        """
        Write out pipe data.
        """
        for (u, v) in self._nw_container.network.edges():
            edge_data = self._nw_container.network[u][v][0]["data"]
            if isinstance(edge_data, PipeContainer):
                dictionary = edge_data.__dict__
                pipe_id = dictionary.get("arc_id")
                if edge_data.flow_fixed:
                    self._coverage_picture(dictionary, pipe_id, u, v)
                for entry, value in dictionary.items():
                    if isinstance(value, unit.Quantity):
                        dictionary[entry] = value.to_tuple()
                    if entry in ["sampling_points_pending", "sampling_points_build_in", "successful_points"]:
                        sps = []
                        for sp in value:
                            sp_dictionary = sp.__dict__
                            sp_dictionary.update({"id": id(sp)})
                            sps.append(sp_dictionary)
                        dictionary[entry] = sps
                with open(str(self._pipe_dir.joinpath(f"{pipe_id}.json")), "w") as f:
                    f.write(json.dumps(dictionary, indent=4))

    def _coverage_picture(self, dictionary: dict, pipe_id: str, u: str, v: str):
        """
        Make a nice coverage picture of the pipes for the 2D-case, that is, for the case where the flow is fixed for
        the pipe.
        """
        u_data = self._nw_container.network.nodes[u]["data"]
        v_data = self._nw_container.network.nodes[v]["data"]
        p_in_min = u_data.pressure_min.to(pressure_unit).magnitude
        p_in_max = u_data.pressure_max.to(pressure_unit).magnitude
        p_out_min = v_data.pressure_min.to(pressure_unit).magnitude
        p_out_max = v_data.pressure_max.to(pressure_unit).magnitude
        fig, ax = plt.subplots()
        p_in_values = np.linspace(p_in_min, p_in_max, 100)
        # We can use flow min, as the flow is fixed
        flow = dictionary["flow_min"].to(flow_unit).magnitude
        mass_flow = norm_vol_per_hour_to_mass_flow(flow)
        coefficient = dictionary["coefficient"]
        p_out_values = np.sqrt(p_in_values ** 2 - coefficient * abs(mass_flow) * mass_flow)
        ax.set_xlim(p_in_min * 0.9, p_in_max * 1.1)
        ax.set_ylim(p_out_min * 0.9, p_out_max * 1.1)
        # ax.xaxis.set_ticks(np.arange(40, 80, 10))
        # ax.yaxis.set_ticks(np.arange(40, 80, 10))
        ax.set_aspect("equal")
        ax.add_patch(
            Rectangle((p_in_min, p_out_min), (p_in_max - p_in_min), (p_out_max - p_out_min), color="royalblue"),
        )
        plt.xlabel(r"$p_u$")
        plt.ylabel(r"$p_v$")
        plt.title(f"{pipe_id}")
        plt.plot(p_in_values, p_out_values, "salmon")
        for entry, value in dictionary.items():
            if entry in ["sampling_points_build_in", "successful_points"]:
                for i, sp in enumerate(value):
                    sp_dictionary = sp.__dict__
                    p_in = sp_dictionary["pressure_in"]
                    p_out = sp_dictionary["pressure_out"]
                    p_in_lb = sp_dictionary["pressure_in_lb"]
                    p_in_ub = sp_dictionary["pressure_in_ub"]
                    p_out_lb = sp_dictionary["pressure_out_lb"]
                    p_out_ub = sp_dictionary["pressure_out_ub"]
                    plt.plot(p_in, p_out, "ok", linewidth=2, markersize=3)
                    # ax.annotate(str(sp_dictionary["iteration"]), (p_in, p_out))
                    ax.add_patch(
                        Rectangle(
                            (p_in_lb, p_out_lb),
                            (p_in_ub - p_in_lb),
                            (p_out_ub - p_out_lb),
                            color="gold",
                            alpha=0.4,
                        ),
                    )
                    # fig.savefig(str(self._pipe_dir.joinpath(f"{pipe_id}_{i}.pdf")), bbox_inches='tight')
        fig.savefig(str(self._pipe_dir.joinpath(f"{pipe_id}.pdf")), bbox_inches="tight")

    def _write_box_solution(self):
        """
        Write out solution.
        """
        if self._model.Status == GRB.OPTIMAL:
            node_to_pressure = {}
            arcs_to_flows = {}
            bnd_nodes_to_flows = {}
            compressors_to_press_inc = {}
            for key, value in self._node_press_vars.items():
                node_to_pressure.update({key: unit.Quantity(value.x, "bar")})
            for key, value in self._edge_flow_vars.items():
                norm_vol = mass_flow_to_norm_vol_per_hour(value.x)
                arcs_to_flows.update({key: unit.Quantity(norm_vol, flow_unit)})
            for key, value in self._node_flow_vars.items():
                norm_vol = mass_flow_to_norm_vol_per_hour(value.x)
                bnd_nodes_to_flows.update({key: unit.Quantity(norm_vol, flow_unit)})
            for key, value in self._compr_vars.items():
                compressors_to_press_inc.update({key: unit.Quantity(value.x, "bar")})
            self._write_solution(node_to_pressure, arcs_to_flows, bnd_nodes_to_flows, compressors_to_press_inc)

    def _write_config(self, time_limit: bool):
        """
        Write out config.
        """
        step_pressure = self._sampling_steps[0] if self._sampling_steps else None
        step_flow = self._sampling_steps[1] if self._sampling_steps else None
        dic = {
            "error_value": self._error_value,
            "absolute_error": self._absolute_error,
            "sampling_steps_pressure_bar": step_pressure,
            "sampling_steps_mass_flow_kg_per_second": step_flow,
            "time_limit": time_limit,
        }
        with open(str(self._result_directory.joinpath("config.json")), "w") as f:
            f.write(json.dumps(dic, indent=4))

    def _write_summary(
        self,
        iteration: int,
        obj: int,
        max_error: float,
        algorithm_runtime: float,
        solver_runtime: float,
    ):
        """
        Write out summary.
        """
        dic = {
            "read_network_runtime": self._read_network_runtime,
            "sampling_runtime_inside_of_read_network_runtime": self._nw_container.sampling_runtime,
            "algorithm_runtime": algorithm_runtime,
            "solver_runtime_inside_of_algorithm_runtime": solver_runtime,
            "total_runtime": self._read_network_runtime + algorithm_runtime,
            "iterations": iteration,
            "max_error": max_error,
            "objective": obj,
        }
        with open(str(self._result_directory.joinpath("summary.json")), "w") as f:
            f.write(json.dumps(dic, indent=4))

    def _write_solution_stats(self, iteration: int) -> int:
        """
        Write out current solution stats and return the runtime in seconds.
        """
        if self._model.Status == GRB.OPTIMAL:
            self._model.write(str(self._sol_dir.joinpath(f"sol_{iteration}.sol")))
        if self._model.Status == GRB.INFEASIBLE:
            self._model.computeIIS()
            self._model.write(str(self._sol_dir.joinpath(f"ilp_{iteration}.ilp")))
        with open(str(self._sol_dir.joinpath(f"solution_stats_{iteration}.json")), "w") as f:
            opt = "optimal"
            inf = "infeasible"
            unkn = "unknown"
            text = opt if self._model.Status == GRB.OPTIMAL else inf if self._model.Status == GRB.INFEASIBLE else unkn
            runtime_in_seconds = self._model.getAttr(GRB.Attr.Runtime)
            f.write(
                json.dumps(
                    {
                        "Status": self._model.Status,
                        "StatusText": text,
                        "Runtime": runtime_in_seconds,
                        "NumVars": self._model.getAttr(GRB.Attr.NumVars),
                        "NumConstrs": self._model.getAttr(GRB.Attr.NumConstrs),
                        "NumBinVars": self._model.getAttr(GRB.Attr.NumBinVars),
                        "NodeCount": self._model.getAttr(GRB.Attr.NodeCount),
                    },
                    indent=4,
                ),
            )
            return runtime_in_seconds
